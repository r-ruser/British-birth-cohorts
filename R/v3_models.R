options(stringsAsFactors=FALSE);library(survey);library(sandwich);library(geepack)
options(survey.lonely.psu='adjust')
old<-'E:/cohort_rds_stage_20260925/life_course_v2_20260928';out<-'E:/cohort_rds_stage_20260925/life_course_v3_20260928';dir.create(out,recursive=TRUE,showWarnings=FALSE)
common<-c('male','low_birthweight','ses','pregnancy_continued_smoke');mcsvars<-c(common,'respondent_asthma','atopy3','bmi3_z','ever_asthma3');adultvars<-c(common,'childhood_symptom')
num<-function(x)suppressWarnings(as.numeric(unclass(x)))
data<-setNames(lapply(c('MCS','NCDS','BCS70'),function(c)readRDS(file.path(old,paste0(c,'_v2_dataset.rds')))),c('MCS','NCDS','BCS70'))
models<-list();modeldiag<-list();flows<-list()
estimate<-function(z,form,c){
 if(c=='MCS'){ds<-svydesign(ids=~family,strata=~design_stratum,weights=~w,data=z,nest=TRUE);f<-svyglm(form,ds,family=quasipoisson('log'));V<-vcov(f)}else{
 z<-droplevels(z[order(z$id,if('age'%in%names(z))z$age else seq_len(nrow(z))),]);z$gid<-match(z$id,unique(z$id));if(any(abs(z$y-round(z$y))>1e-8)){f<-glm(form,data=z,weights=w,family=quasipoisson('log'));V<-vcovCL(f,cluster=z$id,type='HC0')}else{f<-geeglm(form,data=z,id=gid,weights=w,family=poisson('log'),corstr='independence');V<-vcov(f)}
 };list(f=f,V=V)
}
extract<-function(ff,c,endpoint,method,terms=NULL){b<-coef(ff$f);V<-ff$V;if(!is.null(terms)){b<-b[terms];V<-V[terms,terms,drop=FALSE]};s<-sqrt(diag(V));data.frame(cohort=c,endpoint,method,term=names(b),log_RR=b,SE=s,RR=exp(b),lower=exp(b-1.96*s),upper=exp(b+1.96*s),p=2*pnorm(-abs(b/s)))}
collect<-function(ff,c,endpoint,method,terms=NULL){models[[paste(c,endpoint,method,sep='_')]]<<-extract(ff,c,endpoint,method,terms);off<-model.offset(model.frame(ff$f));mu<-fitted(ff$f);if(!is.null(off))mu<-mu/exp(off);modeldiag[[paste(c,endpoint,method,sep='_')]]<<-data.frame(cohort=c,endpoint,method,minimum_fitted=min(mu),maximum_fitted=max(mu),above_one=sum(mu>1))}
ipcw_long<-function(d,ages,vars,c,known=FALSE){
 req<-vars;if(c=='MCS')req<-c(req,'basewt','family','design_stratum');d<-d[complete.cases(d[req]),];if(known)d<-d[!is.na(d[[paste0('y',ages[1])]]),]
 d$base<-if(c=='MCS')d$basewt else 1;if('child_selection_weight'%in%names(d))d$base<-d$base*d$child_selection_weight
 active<-seq_len(nrow(d));cw<-rep(1,nrow(d));last<-hist<-rep(0,nrow(d));zz<-list()
 for(j in seq_along(ages)){age<-ages[j];dead<-if(c=='NCDS')d[[paste0('dead',age)]]else rep(FALSE,nrow(d));active<-active[!dead[active]];y<-d[[paste0('y',age)]];r<-as.integer(!is.na(y[active]));q<-d[active,];q$r<-r;q$last<-last[active];q$history<-hist[active]/max(j-1,1)
  if(known&&j==1){pd<-pn<-rep(1,length(r))}else{fd<-glm(reformulate(c(vars,if(j>1)c('last','history')),'r'),q,family=binomial());pd<-pmin(.98,pmax(.02,predict(fd,type='response')));pn<-rep(mean(r),length(r))}
  cw[active]<-cw[active]*pn/pd;keep<-active[r==1];z<-d[keep,];z$age<-age;z$y<-y[keep];z$previous<-last[keep];z$ipcw_raw<-cw[keep];cuts<-quantile(z$ipcw_raw,c(.01,.99));z$ipcw<-pmin(cuts[2],pmax(cuts[1],z$ipcw_raw));z$w<-z$base*z$ipcw;zz[[j]]<-z
  last[keep]<-y[keep];hist[keep]<-hist[keep]+y[keep];active<-keep
 }
 z<-do.call(rbind,zz);z$wave<-factor(z$age,levels=ages);all<-do.call(rbind,lapply(ages,function(age){u<-d[!is.na(d[[paste0('y',age)]]),];u$age<-age;u$y<-u[[paste0('y',age)]];u$w<-u$base;u$wave<-factor(age,levels=ages);u}));list(monotone=z,all=all,baseline=d)
}
if(Sys.getenv('V3_HELPERS_ONLY')=='1')quit(save='no')
# 1. Equivalence on log-RR differences, independent cohorts.
r<-read.csv(file.path(old,'model_RR_results.csv'));r<-r[r$endpoint=='shared_age_10_11'&r$method=='response_IPW'&r$term%in%common,];eq<-list();kk<-0
for(v in common)for(pair in list(c('MCS','NCDS'),c('MCS','BCS70'),c('NCDS','BCS70'))){a<-r[r$cohort==pair[1]&r$term==v,];b<-r[r$cohort==pair[2]&r$term==v,];delta<-log(a$estimate/b$estimate);se<-sqrt(((log(a$upper)-log(a$lower))/3.92)^2+((log(b$upper)-log(b$lower))/3.92)^2)
 for(bound in c('primary','strict')){lims<-if(bound=='primary')c(.80,1.25)else c(.85,1.18);pL<-pnorm(-(delta-log(lims[1]))/se);pU<-pnorm((delta-log(lims[2]))/se);kk<-kk+1;eq[[kk]]<-data.frame(term=v,comparison=paste(pair,collapse='/'),bound,ROR=exp(delta),lower90=exp(delta-qnorm(.95)*se),upper90=exp(delta+qnorm(.95)*se),margin_lower=lims[1],margin_upper=lims[2],p_lower=pL,p_upper=pU,p_TOST=max(pL,pU),equivalent=max(pL,pU)<.05)}
}
eq<-do.call(rbind,eq);write.csv(eq,file.path(out,'cohort_equivalence_pairwise.csv'),row.names=FALSE)
eg<-aggregate(p_TOST~term+bound,eq,max);eg$all_three_equivalent<-eg$p_TOST<.05;eg$p_Holm_four<-ave(eg$p_TOST,eg$bound,FUN=function(x)p.adjust(x,'holm'));write.csv(eg,file.path(out,'cohort_equivalence_global.csv'),row.names=FALSE)
# 2. Stage differences: fully interacted wave supports direct comparison with separate fits.
ll<-readRDS(file.path(old,'analysis_longitudinal_frames.rds'));z<-droplevels(ll$MCS$monotone[ll$MCS$monotone$age>3,]);st<-list();joint<-list()
for(method in c('stage_specific_wave','shared_wave')){form<-if(method=='stage_specific_wave')as.formula(paste('y ~ previous * (wave +',paste(mcsvars,collapse='+'),')'))else as.formula(paste('y ~ wave + previous * (',paste(mcsvars,collapse='+'),')'))
 ff<-estimate(z,form,'MCS');collect(ff,'MCS','transition_stage',method)
 for(v in mcsvars){nm<-intersect(c(paste0('previous:',v),paste0(v,':previous')),names(coef(ff$f)));stopifnot(length(nm)==1);delta<--coef(ff$f)[nm];se<-sqrt(ff$V[nm,nm]);st[[paste(method,v)]]<-data.frame(method,term=v,onset_RR=exp(coef(ff$f)[v]),persistence_RR=exp(coef(ff$f)[v]+coef(ff$f)[nm]),stage_ROR=exp(delta),lower=exp(delta-1.96*se),upper=exp(delta+1.96*se),p_interaction=2*pnorm(-abs(delta/se)))}
 jj<-which(names(coef(ff$f))%in%c(paste0('previous:',mcsvars),paste0(mcsvars,':previous')));b<-coef(ff$f)[jj];V<-ff$V[jj,jj];chi<-as.numeric(t(b)%*%solve(V,b));joint[[method]]<-data.frame(method,df=length(jj),chi2=chi,p=pchisq(chi,length(jj),lower.tail=FALSE))
}
st<-do.call(rbind,st);st$q_BH<-ave(st$p_interaction,st$method,FUN=function(x)p.adjust(x,'BH'));write.csv(st,file.path(out,'stage_interactions.csv'),row.names=FALSE);write.csv(do.call(rbind,joint),file.path(out,'stage_joint_tests.csv'),row.names=FALSE)
# 5. Common-age absolute risks, empirical weighted covariate standardization.
ar<-list();commonfits<-list()
for(c in names(data)){d<-data[[c]];age<-if(c=='BCS70')10 else 11;req<-common;if(c=='MCS')req<-c(req,'basewt','family','design_stratum');d<-d[complete.cases(d[req]),];if(c=='NCDS')d<-d[!d$dead11,];d$r<-as.integer(!is.na(d[[paste0('y',age)]]));fr<-glm(reformulate(common,'r'),d,family=binomial());pr<-pmin(.98,pmax(.02,predict(fr,type='response')));d$rw<-mean(d$r)/pr;u<-d[d$r==1,];cut<-quantile(u$rw,c(.01,.99));u$w<-pmin(cut[2],pmax(cut[1],u$rw))*(if(c=='MCS')u$basewt else 1);u$y<-u[[paste0('y',age)]];ff<-estimate(u,reformulate(common,'y'),c);b<-coef(ff$f);V<-ff$V;values<-list()
 for(sex in 0:1){nd<-u;nd$male<-sex;X<-model.matrix(reformulate(common),nd);mu<-as.numeric(exp(X%*%b));pw<-u$w/sum(u$w);risk<-sum(pw*mu);g<-colSums(X*(pw*mu));values[[sex+1]]<-list(risk=risk,g=g);se<-sqrt(as.numeric(t(g)%*%V%*%g));ar[[paste(c,sex)]]<-data.frame(cohort=c,measure=if(sex==0)'female_risk'else'male_risk',estimate=risk,lower=exp(log(risk)-1.96*se/risk),upper=exp(log(risk)+1.96*se/risk),n=nrow(u))}
 gd<-values[[2]]$g-values[[1]]$g;rd<-values[[2]]$risk-values[[1]]$risk;se<-sqrt(as.numeric(t(gd)%*%V%*%gd));ar[[paste(c,'RD')]]<-data.frame(cohort=c,measure='risk_difference',estimate=rd,lower=rd-1.96*se,upper=rd+1.96*se,n=nrow(u));gr<-values[[2]]$g/values[[2]]$risk-values[[1]]$g/values[[1]]$risk;lr<-log(values[[2]]$risk/values[[1]]$risk);se<-sqrt(as.numeric(t(gr)%*%V%*%gr));ar[[paste(c,'RR')]]<-data.frame(cohort=c,measure='risk_ratio',estimate=exp(lr),lower=exp(lr-1.96*se),upper=exp(lr+1.96*se),n=nrow(u));commonfits[[c]]<-list(f=ff,standardization_n=nrow(u))
}
write.csv(do.call(rbind,ar),file.path(out,'common_age_absolute_risks.csv'),row.names=FALSE)
# 3/4. Multi-window childhood phenotypes and adult cumulative burden.
childraw<-readRDS(file.path(out,'child_selected_sources.rds'))
adults<-list();tabulations<-list();burdenflow<-list();gradient<-list()
for(c in c('NCDS','BCS70')){d<-data[[c]];src<-childraw[[if(c=='NCDS')'NCDS:ncds0123.dta'else'BCS70:bcs7016x.dta']];sid<-if(c=='NCDS')src$ncdsid else src$bcsid;ever<-num(src[[if(c=='NCDS')'n2617'else'ob9_1']]);recent<-num(src[[if(c=='NCDS')'n2622'else'ob9_18']]);y16<-rep(NA_real_,nrow(src));y16[which(recent%in%1:6|ever==2)]<-0;y16[which(recent==7)]<-1;d$ch16<-y16[match(d$id,as.character(sid))];d$chfirst<-d[[if(c=='NCDS')'y11'else'y10']];d$ch7ever<-d$childhood_symptom
 for(definition in if(c=='NCDS')c('two_recent_windows','three_mixed_windows')else'two_recent_windows'){
  chvars<-if(definition=='three_mixed_windows')c('ch7ever','chfirst','ch16')else c('chfirst','ch16');b<-d[complete.cases(d[common]),];b$child_complete<-as.integer(complete.cases(b[chvars]));sel<-glm(reformulate(common,'child_complete'),b,family=binomial());ps<-pmin(.98,pmax(.02,predict(sel,type='response')));b$child_selection_weight<-mean(b$child_complete)/ps;b<-b[b$child_complete==1,];cut<-quantile(b$child_selection_weight,c(.01,.99));b$child_selection_weight<-pmin(cut[2],pmax(cut[1],b$child_selection_weight));count<-rowSums(b[chvars]);b$phenotype<-factor(ifelse(count==0,'never',ifelse(count==1,'isolated','repeated')),levels=c('never','isolated','repeated'));b$score<-as.integer(b$phenotype)-1
  ages<-if(c=='NCDS')c(23,33,42,50,55,62)else c(26,30,34,38,42,46);vars<-c(common,'phenotype');l<-ipcw_long(b,ages,vars,c)
  for(method in c('all_available','monotone_unweighted','IPCW')){zz<-if(method=='all_available')l$all else l$monotone;if(method=='monotone_unweighted')zz$w<-zz$base;ff<-estimate(zz,reformulate(c('wave',vars),'y'),c);collect(ff,c,paste0('child_phenotype_adult_',definition),method,c(common,'phenotypeisolated','phenotyperepeated'));bbcoef<-coef(ff$f);ii<-c('phenotypeisolated','phenotyperepeated');delta<-bbcoef[ii[2]]-bbcoef[ii[1]];se<-sqrt(ff$V[ii[1],ii[1]]+ff$V[ii[2],ii[2]]-2*ff$V[ii[1],ii[2]]);gradient[[paste(c,definition,method)]]<-data.frame(cohort=c,definition,method,RR_repeated_vs_isolated=exp(delta),lower=exp(delta-1.96*se),upper=exp(delta+1.96*se),p=2*pnorm(-abs(delta/se)));fftrend<-estimate(zz,reformulate(c('wave',common,'score'),'y'),c);collect(fftrend,c,paste0('phenotype_trend_',definition),method,'score')}
  # Complete six waves is the comparable primary cumulative burden sample.
  last<-l$monotone[l$monotone$age==max(ages),];last$n_observed<-6;last$positive_count<-rowSums(last[paste0('y',ages)]);last$positive_proportion<-last$positive_count/6;last$adult_category<-factor(ifelse(last$positive_count==0,'none',ifelse(last$positive_count==1,'isolated','repeated')),levels=c('none','isolated','repeated'));last$y<-as.integer(last$positive_count>=2)
  for(method in c('IPCW','unweighted_complete_six')){u<-last;if(method=='unweighted_complete_six')u$w<-1;ff<-estimate(u,reformulate(vars,'y'),c);collect(ff,c,paste0('adult_repeated_six_',definition),method,c(common,'phenotypeisolated','phenotyperepeated'));u$y<-u$positive_count;ff<-estimate(u,as.formula(paste('y ~',paste(vars,collapse='+'),'+ offset(log(n_observed))')),c);collect(ff,c,paste0('adult_positive_fraction_six_',definition),method,c('phenotypeisolated','phenotyperepeated'))}
  for(ph in levels(last$phenotype)){u<-last[last$phenotype==ph,];for(ac in levels(last$adult_category)){tabulations[[paste(c,definition,ph,ac)]]<-data.frame(cohort=c,definition,phenotype=ph,adult_category=ac,n=sum(u$adult_category==ac),total=nrow(u),weighted_probability=sum(u$w*(u$adult_category==ac))/sum(u$w),weighted_positive_fraction=weighted.mean(u$positive_proportion,u$w))}}
  pp<-l$baseline;pp$n_observed<-rowSums(!is.na(pp[paste0('y',ages)]));pp$positive_count<-rowSums(pp[paste0('y',ages)],na.rm=TRUE);pp<-pp[pp$n_observed>=3,];pp$y<-pp$positive_count;pp$w<-pp$child_selection_weight;ff<-estimate(pp,as.formula(paste('y ~',paste(vars,collapse='+'),'+ offset(log(n_observed))')),c);collect(ff,c,paste0('adult_positive_fraction_ge3_',definition),'observed_window_sensitivity',c('phenotypeisolated','phenotyperepeated'));pp$y<-pp$positive_count/pp$n_observed;ff<-estimate(pp,reformulate(vars,'y'),c);collect(ff,c,paste0('adult_person_mean_fraction_ge3_',definition),'observed_window_sensitivity',c('phenotypeisolated','phenotyperepeated'))
  burdenflow[[paste(c,definition)]]<-data.frame(cohort=c,definition,baseline_core_complete=sum(complete.cases(d[common])),child_all_windows=nrow(b),adult_any=length(unique(l$all$id)),adult_first=nrow(l$monotone[l$monotone$age==min(ages),]),adult_complete_six=nrow(last),adult_ge3=nrow(pp),six_events=sum(last$positive_count>=2),six_ESS=sum(last$w)^2/sum(last$w^2))
  adults[[paste(c,definition,sep='_')]]<-list(frames=l,cumulative_six=last,observed_ge3=pp)
 }
 data[[c]]<-d
}
write.csv(do.call(rbind,burdenflow),file.path(out,'burden_sample_flow.csv'),row.names=FALSE);write.csv(do.call(rbind,tabulations),file.path(out,'child_to_adult_burden_table.csv'),row.names=FALSE)
write.csv(do.call(rbind,gradient),file.path(out,'child_phenotype_gradient_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,models),file.path(out,'new_model_results.csv'),row.names=FALSE);write.csv(do.call(rbind,modeldiag),file.path(out,'new_model_diagnostics.csv'),row.names=FALSE)
saveRDS(adults,file.path(out,'burden_analysis_frames.rds'));saveRDS(data,file.path(out,'v3_cohort_datasets.rds'));saveRDS(commonfits,file.path(out,'common_age_models.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));cat('V3 principal models complete\n')
