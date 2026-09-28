out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
library(survey);library(sandwich)
options(survey.lonely.psu='adjust')
exposures<-c('pregnancy_continued_smoke','mother_smoke','father_smoke','any_parent_smoke','breastfed')
results<-list();flows<-list();diagnostics<-list();risks<-list();k<-0L;kd<-0L;kr<-0L
record <- function(fit,V,cohort,endpoint,level,exposure,n,events,rows,scale) {
 b<-coef(fit)[exposure];se<-sqrt(diag(V))[exposure]
 data.frame(cohort,endpoint,adjustment=level,exposure,scale,estimate=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p=2*pnorm(-abs(b/se)),people=n,events,interval_rows=rows,log_estimate=b,SE=se)
}
mcperiods <- function(d,base,ages,stem,vars) {
 d<-d[!is.na(d[[base]]) & d[[base]]==0 & complete.cases(d[vars]),]
 zz<-list();kk<-0L;active<-seq_len(nrow(d));previous<-3
 for(age in ages) {
  y<-d[[paste0(stem,age)]];wt<-d[[paste0('weight_age',age)]]
  active<-active[!is.na(y[active]) & !is.na(wt[active]) & !is.na(d$design_stratum[active])]
  if(!length(active)) break
  kk<-kk+1L;zz[[kk]]<-data.frame(id=d$id[active],family=d$family[active],stratum=d$design_stratum[active],event=y[active],
    interval=paste0(previous,'-',age),duration=age-previous,weight=wt[active],d[active,vars,drop=FALSE],check.names=FALSE)
  active<-active[y[active]==0];previous<-age
 }
 z<-do.call(rbind,zz);z$interval<-factor(z$interval,levels=unique(z$interval));z
}
for(cohort in c('MCS','NCDS','BCS70')) {
 d<-readRDS(file.path(out,paste0(cohort,'_parent_behavior_dataset.rds')))
 for(exposure in exposures) {
  if(sum(!is.na(d[[exposure]]))==0) next
  for(level in c('basic','extended')) {
   covars<-c('male','ses')
   if(exposure!='pregnancy_continued_smoke') covars<-c(covars,'pregnancy_continued_smoke')
   if(level=='extended') covars<-c(covars,if(cohort=='MCS') c('low_birthweight','room_smoke','respondent_asthma','atopy3','bmi3_z') else c('low_birthweight','childhood_symptom'))
   vars<-unique(c(covars,exposure))
   if(cohort=='MCS') {
    # Baseline current/ever symptoms added only in extended conditional models.
    endpoints<-c('first_asthma_report','first_wheeze_report','wheeze_age11_after_age3_wheeze')
    for(endpoint in endpoints) {
     vv<-vars
     if(level=='extended') vv<-unique(c(vv,if(endpoint=='first_asthma_report') 'ever_wheeze3' else if(endpoint=='wheeze_age11_after_age3_wheeze') 'ever_asthma3' else character()))
     if(endpoint=='wheeze_age11_after_age3_wheeze') {
      z<-d[!is.na(d$wheeze3) & d$wheeze3==1 & !is.na(d$wheeze11) & complete.cases(d[c(vv,'weight_age11','family','design_stratum')]),]
      z$event<-z$wheeze11;z$weight<-z$weight_age11;z$stratum<-z$design_stratum
      form<-reformulate(vv,'event');fam<-quasipoisson(link='log');scale<-'RR'
     } else {
      if(endpoint=='first_asthma_report') z<-mcperiods(d,'ever_asthma3',c(5,7,11,14),'ever_asthma',vv)
      else z<-mcperiods(d,'ever_wheeze3',c(5,7,11),'wheeze',vv)
      form<-as.formula(paste('event ~ interval +',paste(vv,collapse='+'),'+ offset(log(duration))'))
      fam<-quasibinomial(link='cloglog');scale<-'HR'
     }
     stopifnot(all(z$weight>0),!anyNA(z$event))
     des<-svydesign(ids=~family,strata=~stratum,weights=~weight,data=z,nest=TRUE)
     fit<-svyglm(form,design=des,family=fam)
     k<-k+1L;results[[k]]<-record(fit,vcov(fit),cohort,endpoint,level,exposure,length(unique(z$id)),sum(z$event),nrow(z),scale)
     saveRDS(list(model=fit,exposure=exposure,formula=form),file.path(out,paste(cohort,endpoint,level,exposure,'model.rds',sep='_')),compress='xz')
     if(scale=='HR') {
      alt<-svyglm(update(form,as.formula(paste('. ~ . + interval:',exposure))),design=des,family=fam)
      tst<-regTermTest(alt,as.formula(paste('~ interval:',exposure)))
      kd<-kd+1L;diagnostics[[kd]]<-data.frame(cohort,endpoint,adjustment=level,exposure,time_interaction_p=tst$p,converged=fit$converged)
      if(level=='basic') for(g in 0:1) {
       a<-z[z[[exposure]]==g,];h<-sapply(levels(z$interval),function(iv){b<-a[a$interval==iv,];weighted.mean(b$event,b$weight)})
       kr<-kr+1L;risks[[kr]]<-data.frame(cohort,endpoint,exposure,level=g,people=length(unique(a$id)),events=sum(a$event),crude_weighted_cumulative_risk=1-prod(1-h))
      }
     } else {
      kd<-kd+1L;diagnostics[[kd]]<-data.frame(cohort,endpoint,adjustment=level,exposure,time_interaction_p=NA_real_,converged=fit$converged)
      if(level=='basic') for(g in 0:1) {
       a<-z[z[[exposure]]==g,];kr<-kr+1L;risks[[kr]]<-data.frame(cohort,endpoint,exposure,level=g,people=nrow(a),events=sum(a$event),crude_weighted_cumulative_risk=weighted.mean(a$event,a$weight))
      }
     }
     cat(cohort,endpoint,level,exposure,length(unique(z$id)),sum(z$event),'\n');flush.console()
    }
   } else {
    ages<-if(cohort=='NCDS') c(50,55,62) else c(38,42,46)
    ys<-paste0('y',ages);a<-d[complete.cases(d[vars]),]
    response_model<-function(a,y) {
     a$response<-as.integer(!is.na(a[[y]]))
     f<-glm(reformulate(vars,'response'),data=a,family=binomial())
     pp<-pmin(pmax(predict(f,type='response'),.05),.95)
     mean(a$response)/pp
    }
    a$w0<-response_model(a,ys[1]);b<-a[!is.na(a[[ys[1]]]) & a[[ys[1]]]==0,]
    b$w1<-b$w0*response_model(b,ys[2]);c<-b[!is.na(b[[ys[2]]]) & b[[ys[2]]]==0,]
    c$w2<-c$w1*response_model(c,ys[3])
    r1<-b[!is.na(b[[ys[2]]]),c('id',vars,'w1',ys[2])];names(r1)[names(r1)=='w1']<-'weight';names(r1)[names(r1)==ys[2]]<-'event'
    r1$interval<-paste0(ages[1],'-',ages[2]);r1$duration<-diff(ages)[1]
    r2<-c[!is.na(c[[ys[3]]]),c('id',vars,'w2',ys[3])];names(r2)[names(r2)=='w2']<-'weight';names(r2)[names(r2)==ys[3]]<-'event'
    r2$interval<-paste0(ages[2],'-',ages[3]);r2$duration<-diff(ages)[2]
    z<-rbind(r1,r2);z$interval<-factor(z$interval,levels=unique(z$interval));qq<-quantile(z$weight,c(.01,.99));z$weight<-pmin(pmax(z$weight,qq[1]),qq[2])
    form<-as.formula(paste('event ~ interval +',paste(vars,collapse='+'),'+ offset(log(duration))'))
    fit<-glm(form,data=z,weights=weight,family=quasibinomial(link='cloglog'))
    V<-vcovCL(fit,cluster=z$id,type='HC0');endpoint<-'adult_first_rereport'
    k<-k+1L;results[[k]]<-record(fit,V,cohort,endpoint,level,exposure,length(unique(z$id)),sum(z$event),nrow(z),'HR')
    saveRDS(list(model=fit,robust_vcov=V,exposure=exposure,formula=form,weight_bounds=qq),file.path(out,paste(cohort,endpoint,level,exposure,'model.rds',sep='_')),compress='xz')
    alt<-glm(update(form,as.formula(paste('. ~ . + interval:',exposure))),data=z,weights=weight,family=quasibinomial(link='cloglog'))
    AV<-vcovCL(alt,cluster=z$id,type='HC0');nm<-grep(paste0('interval.*:',exposure,'|',exposure,':interval'),names(coef(alt)),value=TRUE)
    pv<-2*pnorm(-abs(coef(alt)[nm]/sqrt(diag(AV))[nm]))
    kd<-kd+1L;diagnostics[[kd]]<-data.frame(cohort,endpoint,adjustment=level,exposure,time_interaction_p=as.numeric(pv),converged=fit$converged)
    if(level=='basic') for(g in 0:1) {
     a<-z[z[[exposure]]==g,];h<-sapply(levels(z$interval),function(iv){t<-a[a$interval==iv,];weighted.mean(t$event,t$weight)})
     kr<-kr+1L;risks[[kr]]<-data.frame(cohort,endpoint,exposure,level=g,people=length(unique(a$id)),events=sum(a$event),crude_weighted_cumulative_risk=1-prod(1-h))
    }
    cat(cohort,level,exposure,length(unique(z$id)),sum(z$event),'\n');flush.console()
   }
  }
 }
}
r<-do.call(rbind,results);r$q_BH<-NA_real_
for(f in unique(paste(r$cohort,r$endpoint,r$adjustment))) {
 ix<-which(paste(r$cohort,r$endpoint,r$adjustment)==f);r$q_BH[ix]<-p.adjust(r$p[ix],'BH')
}
write.csv(r,file.path(out,'parent_behavior_HR_RR.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'parent_behavior_model_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,risks),file.path(out,'parent_behavior_crude_absolute_risk.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
