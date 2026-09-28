options(stringsAsFactors=FALSE)
library(survey);library(sandwich);library(geepack)
options(survey.lonely.psu='adjust')
out<-'E:/cohort_rds_stage_20260925/life_course_v2_20260928'
raw<-readRDS(file.path(out,'selected_source_tables.rds'))
num<-function(x)suppressWarnings(as.numeric(unclass(x)))
yn<-function(x,yes=1,no=2)ifelse(num(x)==yes,1,ifelse(num(x)==no,0,NA_real_))
tab<-function(c,f)raw[[paste(c,f,sep=':')]]
idcol<-function(d)names(d)[tolower(names(d)) %in% c('ncdsid','bcsid')][1]
route<-function(screen,detail,yes=1,no=2,sno=2){z<-yn(detail,yes,no);z[!is.na(screen)&num(screen)==sno & is.na(z)]<-0;z}
joiny<-function(a,d,y,age){b<-data.frame(id=as.character(d[[idcol(d)]]),y=y);names(b)[2]<-paste0('y',age);stopifnot(!anyDuplicated(b$id));merge(a,b,by='id',all.x=TRUE,sort=FALSE)}
data<-list()
for(c in c('NCDS','BCS70','MCS'))data[[c]]<-readRDS(paste0('E:/cohort_rds_stage_20260925/parent_behavior_20260928/',c,'_parent_behavior_dataset.rds'))
a<-data$NCDS
vital<-readRDS(file.path(out,'NCDS_response_source.rds'));ix<-match(a$id,as.character(vital$ncdsid))
for(pair in list(c(11,'outcme02'),c(23,'outcme04'),c(33,'outcme05'),c(42,'outcme06'),c(50,'outcme08'),c(55,'outcme09'),c(62,'outcme10'))){status<-num(vital[[pair[2]]])[ix];a[[paste0('dead',pair[1])]]<-ifelse(is.na(status),FALSE,status %in% c(8,10))}
d<-tab('NCDS','ncds4.dta');a<-joiny(a,d,route(d$n5770,d$n5771),23)
d<-tab('NCDS','ncds5cmi.dta');a<-joiny(a,d,route(d$n504018,d$n504019),33)
d<-tab('NCDS','ncds6_v2.dta');a<-joiny(a,d,route(d$wheezy,d$wheezy12),42)
d<-tab('NCDS','ncds0123.dta');y<-ifelse(num(d$n1306)%in%1:4,1,ifelse(num(d$n1306)==5|num(d$n1305)==3,0,NA_real_));a<-joiny(a,d,y,11)
data$NCDS<-a
a<-data$BCS70;d<-tab('BCS70','bcs96x.dta')
# Checkbox: an unselected 12-month item among valid since-16 respondents means no.
y<-rep(NA_real_,nrow(d));y[which(num(d$b960448)==0 | (num(d$b960448)==1 & is.na(d$b960514)))]<-0;y[which(num(d$b960514)==1)]<-1;y[which(num(d$b960514)<0)]<-NA;a<-joiny(a,d,y,26)
d<-tab('BCS70','bcs2000.dta');a<-joiny(a,d,route(d$wheezy,d$wheezy12),30)
d<-tab('BCS70','bcs_2004_followup.dta');a<-joiny(a,d,route(d$bd7hpb01,d$b7asth2m,sno=0),34)
d<-tab('BCS70','sn3723.dta');a<-joiny(a,d,route(d$b7_1,d$b7_20),10);data$BCS70<-a
a<-data$MCS;d<-tab('MCS','mcs_longitudinal_cm_file.dta');key<-paste(as.character(d$mcsid),num(d$cnum),sep=':');a$basewt<-num(d$bovwt2)[match(a$id,key)];a$basewt[a$basewt<=0]<-NA
for(age in c(3,5,7,11))a[[paste0('y',age)]]<-a[[paste0('wheeze',age)]]
data$MCS<-a
for(c in names(data))saveRDS(data[[c]],file.path(out,paste0(c,'_v2_dataset.rds')))
common<-c('male','low_birthweight','ses','pregnancy_continued_smoke')
adultvars<-c(common,'childhood_symptom')
mcsvars<-c(common,'respondent_asthma','atopy3','bmi3_z','ever_asthma3')
results<-list();flows<-list();diag<-list();rates<-list();ints<-list();sexcontrasts<-list();modeldiag<-list();k<-0;ki<-0
record<-function(f,V,c,endpoint,method,n,rows,events){b<-coef(f);s<-sqrt(diag(V));z<-data.frame(cohort=c,endpoint,method,term=names(b),estimate=exp(b),lower=exp(b-1.96*s),upper=exp(b+1.96*s),p=2*pnorm(-abs(b/s)),people=n,rows,events);k<<-k+1;results[[k]]<<-z}
fit<-function(z,formula,c,endpoint,method){stopifnot(all(z$w>0),!anyNA(z$y));if(c=='MCS'){
 ds<-svydesign(ids=~family,strata=~design_stratum,weights=~w,data=z,nest=TRUE);f<-svyglm(formula,ds,family=quasipoisson('log'));V<-vcov(f)
 }else if(endpoint!='shared_age_10_11'){z<-droplevels(z[order(z$id,z$age),]);z$gid<-match(z$id,unique(z$id));f<-geeglm(formula,data=z,family=poisson('log'),weights=w,id=gid,corstr='independence',std.err='san.se');V<-vcov(f)
 }else{f<-glm(formula,z,family=poisson('log'),weights=w);V<-vcovCL(f,cluster=z$id,type='HC0')}
 record(f,V,c,endpoint,method,length(unique(z$id)),nrow(z),sum(z$y));modeldiag[[paste(c,endpoint,method,sep='_')]]<<-data.frame(cohort=c,endpoint,method,min_fitted=min(fitted(f)),max_fitted=max(fitted(f)),fitted_above_one=sum(fitted(f)>1));list(f=f,V=V)
}
joint<-function(ff,c,end,method){jj<-grep(':',names(coef(ff$f)));b<-coef(ff$f)[jj];V<-ff$V[jj,jj,drop=FALSE];chi<-as.numeric(t(b)%*%solve(V,b));ki<<-ki+1;ints[[ki]]<<-data.frame(cohort=c,endpoint=end,method,df=length(jj),chi2=chi,p=pchisq(chi,length(jj),lower.tail=FALSE))}
# Sequential stabilized IPCW: common baseline complete cases; artificial monotone censoring.
longify<-function(d,ages,vars,c,baseline_known=FALSE){
 req<-vars;if(c=='MCS')req<-c(req,'basewt','family','design_stratum')
 d<-d[complete.cases(d[req]),];if(baseline_known)d<-d[!is.na(d[[paste0('y',ages[1])]]),]
 d$bw<-if(c=='MCS')d$basewt else 1;active<-seq_len(nrow(d));cumw<-rep(1,nrow(d));last<-rep(0,nrow(d));hist<-rep(0,nrow(d));zz<-list()
 for(j in seq_along(ages)){
  age<-ages[j];yy<-d[[paste0('y',age)]];dead<-if(c=='NCDS')d[[paste0('dead',age)]] else rep(FALSE,nrow(d));active<-active[!dead[active]];r<-as.integer(!is.na(yy[active]));q<-d[active,];q$r<-r;q$last<-last[active];q$history<-hist[active]/max(j-1,1)
  if(baseline_known && j==1){pd<-pn<-rep(1,length(active))}else{
   denvars<-c(vars,if(j>1)c('last','history'))
   fd<-glm(reformulate(denvars,'r'),q,family=binomial());pd0<-predict(fd,type='response');pd<-pmin(.98,pmax(.02,pd0));pn<-rep(mean(r),length(r))
   diag[[paste(c,age,sep='_')]]<<-data.frame(cohort=c,age,riskset=length(active),observed=sum(r),min_probability=min(pd0),max_probability=max(pd0),clipped_probabilities=sum(pd!=pd0))
  }
  cumw[active]<-cumw[active]*pn/pd
  keep<-active[r==1];z<-d[keep,];z$age<-age;z$y<-yy[keep];z$previous<-last[keep];z$history_previous<-hist[keep];z$ipcw_raw<-cumw[keep];z$baseline_w<-d$bw[keep]
  if(nrow(z)){cut<-quantile(z$ipcw_raw,c(.01,.99));z$ipcw<-pmin(cut[2],pmax(cut[1],z$ipcw_raw));z$w<-z$baseline_w*z$ipcw;zz[[j]]<-z}
  flows[[paste(c,age,sep='_')]]<<-data.frame(cohort=c,age,baseline_complete=nrow(d),known_dead=if(c=='NCDS')sum(dead)else NA_integer_,monotone_observed=length(keep),all_available=sum(!is.na(yy)),events=sum(yy[keep]),ipcw_min=if(length(keep))min(z$ipcw_raw)else NA,ipcw_max=if(length(keep))max(z$ipcw_raw)else NA,ESS=if(length(keep))sum(z$w)^2/sum(z$w^2)else NA)
  last[keep]<-yy[keep];hist[keep]<-hist[keep]+yy[keep];active<-keep
 }
 z<-do.call(rbind,zz);z$wave<-factor(z$age,levels=ages)
 all<-do.call(rbind,lapply(ages,function(age){u<-d[!is.na(d[[paste0('y',age)]]),];u$y<-u[[paste0('y',age)]];u$age<-age;u$w<-u$bw;u$wave<-factor(age,levels=ages);u}))
 list(monotone=z,all=all,baseline=d)
}
longs<-list()
for(c in c('MCS','NCDS','BCS70')){
 ages<-switch(c,MCS=c(3,5,7,11),NCDS=c(23,33,42,50,55,62),BCS70=c(26,30,34,38,42,46));vars<-if(c=='MCS')mcsvars else adultvars
 ll<-longify(data[[c]],ages,vars,c,c=='MCS');longs[[c]]<-ll
 for(method in c('all_available','monotone_unweighted','IPCW','IPCW_untrimmed')){
  z<-if(method=='all_available')ll$all else ll$monotone
  if(method=='monotone_unweighted')z$w<-z$baseline_w
  if(method=='IPCW_untrimmed')z$w<-z$baseline_w*z$ipcw_raw
  endpoint<-if(c=='MCS')'childhood_repeated_wheeze' else 'adult_repeated_burden'
  f<-fit(z,reformulate(c('wave',vars),'y'),c,endpoint,method)
  fi<-fit(z,as.formula(paste('y ~ wave * male +',paste(setdiff(vars,'male'),collapse='+'))),c,'sex_age',method);joint(fi,c,'sex_age',method)
  for(age in ages){contrast<-rep(0,length(coef(fi$f)));names(contrast)<-names(coef(fi$f));contrast['male']<-1;nm<-paste0('wave',age,':male');if(nm %in% names(contrast))contrast[nm]<-1;b<-sum(contrast*coef(fi$f));se<-sqrt(as.numeric(t(contrast)%*%fi$V%*%contrast));sexcontrasts[[paste(c,method,age,sep='_')]]<-data.frame(cohort=c,method,age,RR=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se))}
  for(age in ages)for(sex in 0:1){u<-z[z$age==age & z$male==sex,];if(nrow(u))rates[[paste(c,method,age,sex,sep='_')]]<-data.frame(cohort=c,method,age,male=sex,n=nrow(u),events=sum(u$y),risk=weighted.mean(u$y,u$w))}
  if(c!='MCS'){
   for(window in c('early','late')){u<-z[z$age %in% if(window=='early')ages[1:3] else ages[4:6],];fit(u,reformulate(c('wave',vars),'y'),c,paste0('adult_',window,'_burden'),method)}
  }
 }
}
ll<-longs$MCS
for(method in c('all_adjacent','monotone_unweighted','IPCW','IPCW_untrimmed')){
 if(method=='all_adjacent'){
  d<-ll$baseline;z<-do.call(rbind,lapply(2:4,function(j){ages<-c(3,5,7,11);u<-d[!is.na(d[[paste0('y',ages[j-1])]]) & !is.na(d[[paste0('y',ages[j])]]),];u$previous<-u[[paste0('y',ages[j-1])]];u$y<-u[[paste0('y',ages[j])]];u$wave<-factor(ages[j],levels=c(5,7,11));u$w<-u$bw;u}))
 }else{z<-ll$monotone;z<-z[z$age>3,];z$wave<-droplevels(z$wave);if(method=='monotone_unweighted')z$w<-z$baseline_w;if(method=='IPCW_untrimmed')z$w<-z$baseline_w*z$ipcw_raw}
 for(origin in 0:1){u<-z[z$previous==origin,];fit(u,reformulate(c('wave',mcsvars),'y'),'MCS',if(origin==0)'transition_0_to_1' else 'transition_1_to_1',method)}
 if(method=='IPCW')for(state in c('first_observed','recurrent','symptom_absence')){u<-switch(state,first_observed=z[z$previous==0&z$history_previous==0,],recurrent=z[z$previous==0&z$history_previous>0,],symptom_absence=z[z$previous==1,]);if(state=='symptom_absence')u$y<-1-u$y;fit(u,reformulate(c('wave',mcsvars),'y'),'MCS',paste0('transition_',state),method)}
}
# Shared-age response IPW, fitting identical early-life predictors in each cohort.
replication<-list()
for(c in c('MCS','NCDS','BCS70')){
 age<-if(c=='BCS70')10 else 11;d<-data[[c]];req<-common;if(c=='MCS')req<-c(req,'basewt','family','design_stratum');d<-d[complete.cases(d[req]),];if(c=='NCDS')d<-d[!d$dead11,];d$r<-as.integer(!is.na(d[[paste0('y',age)]]));fr<-glm(reformulate(common,'r'),d,family=binomial());pd<-pmin(.98,pmax(.02,predict(fr,type='response')));d$rw<-mean(d$r)/pd;d$y<-d[[paste0('y',age)]];u<-d[d$r==1,];u$bw<-if(c=='MCS')u$basewt else 1
 for(method in c('unweighted','response_IPW')){u$w<-u$bw;if(method=='response_IPW'){cuts<-quantile(u$rw,c(.01,.99));u$w<-u$bw*pmin(cuts[2],pmax(cuts[1],u$rw))};fit(u,reformulate(common,'y'),c,'shared_age_10_11',method)}
 u$cohort<-c;u$w<-u$w/mean(u$w);u$cluster<-paste(c,if(c=='MCS')u$family else u$id,sep=':');replication[[c]]<-u[c('y','cohort','cluster','w',common)]
}
p<-do.call(rbind,replication);p$cohort<-factor(p$cohort)
for(v in common){f<-glm(as.formula(paste('y ~ cohort *',v,'+',paste(setdiff(common,v),collapse='+'))),p,family=poisson('log'),weights=w);V<-vcovCL(f,cluster=p$cluster,type='HC0');joint(list(f=f,V=V),'three_cohorts',paste0('heterogeneity_',v),'response_IPW')}
res<-do.call(rbind,results);res$q<-NA_real_
sel<-!grepl('wave|Intercept|:',res$term);res$q[sel]<-ave(res$p[sel],interaction(res$endpoint[sel],res$method[sel]),FUN=function(x)p.adjust(x,'BH'))
write.csv(res,file.path(out,'model_RR_results.csv'),row.names=FALSE)
intdf<-do.call(rbind,ints);intdf$q<-ave(intdf$p,interaction(intdf$method,grepl('heterogeneity',intdf$endpoint)),FUN=function(x)p.adjust(x,'BH'));write.csv(intdf,file.path(out,'interaction_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,flows),file.path(out,'attrition_flow_weights.csv'),row.names=FALSE)
write.csv(do.call(rbind,diag),file.path(out,'response_model_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,rates),file.path(out,'sex_age_observed_risks.csv'),row.names=FALSE)
write.csv(do.call(rbind,sexcontrasts),file.path(out,'sex_age_adjusted_RR.csv'),row.names=FALSE)
write.csv(do.call(rbind,modeldiag),file.path(out,'model_fit_diagnostics.csv'),row.names=FALSE)
saveRDS(longs,file.path(out,'analysis_longitudinal_frames.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Complete:',nrow(res),'estimates;',length(ints),'joint tests\n')
