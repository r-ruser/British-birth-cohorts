out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
library(survey);options(survey.lonely.psu='adjust')
ee<-parse('parent_behavior_models.R')
for(e in ee) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e)
d<-readRDS(file.path(out,'MCS_parent_behavior_dataset.rds'))
checks<-read.csv(file.path(out,'parent_behavior_model_checks.csv'))
rr<-read.csv(file.path(out,'parent_behavior_HR_RR.csv'))
stopifnot(nrow(rr)==48,all(checks$converged),all(is.finite(rr$estimate)),all(rr$CI_low>0),all(rr$CI_high>=rr$CI_low))
audit<-list()
for(cc in c('MCS','NCDS','BCS70')) {
 a<-readRDS(file.path(out,paste0(cc,'_parent_behavior_dataset.rds')))
 vv<-c('pregnancy_continued_smoke','mother_smoke','father_smoke','any_parent_smoke','breastfed')
 stopifnot(!anyDuplicated(a$id),all(vapply(a[vv],function(z)all(z[!is.na(z)]%in%0:1),logical(1))))
 audit[[cc]]<-data.frame(cohort=cc,rows=nrow(a),unique_ids=length(unique(a$id)),all_behaviors_binary=TRUE)
}
write.csv(do.call(rbind,audit),file.path(out,'parent_behavior_dataset_checks.csv'),row.names=FALSE)
sel<-checks[checks$cohort=='MCS' & checks$endpoint=='first_asthma_report' & !is.na(checks$time_interaction_p) & checks$time_interaction_p<.05,]
rows<-list();k<-0L
for(i in seq_len(nrow(sel))) {
 s<-sel[i,];obj<-readRDS(file.path(out,paste(s$cohort,s$endpoint,s$adjustment,s$exposure,'model.rds',sep='_')))
 vars<-setdiff(all.vars(obj$formula),c('event','interval','duration'))
 z<-mcperiods(d,'ever_asthma3',c(5,7,11,14),'ever_asthma',vars)
 des<-svydesign(ids=~family,strata=~stratum,weights=~weight,data=z,nest=TRUE)
 fit<-svyglm(update(obj$formula,as.formula(paste('. ~ . + interval:',s$exposure))),design=des,family=quasibinomial(link='cloglog'))
 co<-coef(fit);V<-vcov(fit)
 for(iv in levels(z$interval)) {
  L<-setNames(rep(0,length(co)),names(co));L[s$exposure]<-1
  nn<-intersect(c(paste0('interval',iv,':',s$exposure),paste0(s$exposure,':interval',iv)),names(co))
  if(length(nn))L[nn]<-1
  b<-sum(L*co);se<-sqrt(as.numeric(t(L)%*%V%*%L));k<-k+1L
  rows[[k]]<-data.frame(cohort=s$cohort,endpoint=s$endpoint,adjustment=s$adjustment,exposure=s$exposure,interval=iv,HR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se))
 }
}
write.csv(do.call(rbind,rows),file.path(out,'parent_behavior_interval_specific_HR.csv'),row.names=FALSE)
capture.output(cat('Main estimate rows:',nrow(rr),'\nAll models converged:',all(checks$converged),'\nFinite estimates:',all(is.finite(rr$estimate)),'\n'),file=file.path(out,'analysis_verification.txt'))
