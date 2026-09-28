out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
library(sandwich)
vars<-c('male','maternal_smoke','low_birthweight','low_social','childhood_symptom')
allres<-list();allflow<-list();allrisk<-list();phres<-list()

for(cohort in c('NCDS','BCS70')) {
 d<-readRDS(file.path(out,paste0(cohort,'_adult_analysis_dataset.rds')))
 ages<-if(cohort=='NCDS') c(50,55,62) else c(38,42,46)
 y0<-paste0('y',ages[1]);y1<-paste0('y',ages[2]);y2<-paste0('y',ages[3])
 flow<-data.frame(cohort,stage=c('birth_records','complete_early_predictors','adult_baseline_observed','adult_baseline_no','followup_1_observed_among_baseline_no','followup_2_observed_after_first_no'),
 n=c(nrow(d),sum(complete.cases(d[vars])),sum(!is.na(d[[y0]])),sum(d[[y0]]==0,na.rm=TRUE),sum(d[[y0]]==0 & !is.na(d[[y1]]),na.rm=TRUE),sum(d[[y0]]==0 & d[[y1]]==0 & !is.na(d[[y2]]),na.rm=TRUE)))
 allflow[[cohort]]<-flow
 miss<-data.frame(cohort,variable=c(vars,y0,y1,y2),n_total=nrow(d),n_missing=vapply(d[c(vars,y0,y1,y2)],function(z) sum(is.na(z)),0L))
 miss$missing_percent<-round(100*miss$n_missing/miss$n_total,2)
 write.csv(miss,file.path(out,paste0(cohort,'_adult_missingness.csv')),row.names=FALSE)

 d<-d[complete.cases(d[vars]),]
 # Stabilized response probabilities at adult baseline and at the two later sweeps.
 f0<-glm(as.integer(!is.na(d[[y0]]))~male+maternal_smoke+low_birthweight+low_social+childhood_symptom,data=d,family=binomial())
 p0<-pmin(pmax(predict(f0,type='response'),.05),.95)
 w0<-mean(!is.na(d[[y0]]))/p0
 b<-d[d[[y0]]==0 & !is.na(d[[y0]]),]
 b$w0<-w0[d[[y0]]==0 & !is.na(d[[y0]])]
 f1<-glm(as.integer(!is.na(b[[y1]]))~male+maternal_smoke+low_birthweight+low_social+childhood_symptom,data=b,family=binomial())
 p1<-pmin(pmax(predict(f1,type='response'),.05),.95)
 b$w1<-b$w0*mean(!is.na(b[[y1]]))/p1
 risk2<-b[!is.na(b[[y1]]) & b[[y1]]==0,]
 f2<-glm(as.integer(!is.na(risk2[[y2]]))~male+maternal_smoke+low_birthweight+low_social+childhood_symptom,data=risk2,family=binomial())
 p2<-pmin(pmax(predict(f2,type='response'),.05),.95)
 risk2$w2<-risk2$w1*mean(!is.na(risk2[[y2]]))/p2
 r1<-b[!is.na(b[[y1]]),c('id',vars,'w1',y1)]
 names(r1)[names(r1)=='w1']<-'ipw';names(r1)[names(r1)==y1]<-'event';r1$interval<-paste0(ages[1],'-',ages[2]);r1$duration<-ages[2]-ages[1]
 r2<-risk2[!is.na(risk2[[y2]]),c('id',vars,'w2',y2)]
 names(r2)[names(r2)=='w2']<-'ipw';names(r2)[names(r2)==y2]<-'event';r2$interval<-paste0(ages[2],'-',ages[3]);r2$duration<-ages[3]-ages[2]
 z<-rbind(r1,r2);z$interval<-factor(z$interval,levels=c(unique(r1$interval),unique(r2$interval)))
 q<-quantile(z$ipw,c(.01,.99),na.rm=TRUE);z$ipw_trim<-pmin(pmax(z$ipw,q[1]),q[2])
 saveRDS(z,file.path(out,paste0(cohort,'_adult_risk_intervals.rds')),compress='xz')
 diag<-data.frame(cohort,n_interval_rows=nrow(z),people=length(unique(z$id)),events=sum(z$event),ipw_min=min(z$ipw),ipw_p01=q[1],ipw_median=median(z$ipw),ipw_p99=q[2],ipw_max=max(z$ipw))
 write.csv(diag,file.path(out,paste0(cohort,'_adult_weight_diagnostics.csv')),row.names=FALSE)
 form<-event~interval+male+maternal_smoke+low_birthweight+low_social+childhood_symptom+offset(log(duration))
 for(method in c('unweighted','response_ipw')) {
  fit<-glm(form,data=z,family=quasibinomial(link='cloglog'),weights=if(method=='unweighted') rep(1,nrow(z)) else z$ipw_trim)
  V<-sandwich::vcovCL(fit,cluster=z$id,type='HC0');co<-coef(fit);se<-sqrt(diag(V));j<-names(co)%in%vars
  allres[[paste(cohort,method)]]<-data.frame(cohort,method,predictor=names(co)[j],HR=exp(co[j]),CI_low=exp(co[j]-1.96*se[j]),CI_high=exp(co[j]+1.96*se[j]),p=2*pnorm(-abs(co[j]/se[j])),people=length(unique(z$id)),events=sum(z$event))
  saveRDS(fit,file.path(out,paste0(cohort,'_adult_HR_',method,'.rds')),compress='xz')
 }
 # Check whether the childhood symptom association is stable across the two adult intervals.
 alt<-glm(update(form,. ~ . + interval:childhood_symptom),data=z,family=quasibinomial(link='cloglog'),weights=z$ipw_trim)
 V<-sandwich::vcovCL(alt,cluster=z$id,type='HC0');nm<-grep('interval.*:childhood_symptom|childhood_symptom:interval',names(coef(alt)),value=TRUE)
 v<-coef(alt)[nm];se<-sqrt(diag(V))[nm]
 phres[[cohort]]<-data.frame(cohort,interaction_logHR=v,SE=se,p=2*pnorm(-abs(v/se)))
 # Weighted empirical risk through the final sweep, conditional on an observed negative adult baseline.
 for(g in 0:1) {
  a<-z[z$childhood_symptom==g,]
  h<-sapply(levels(z$interval),function(iv){t<-a[a$interval==iv,];weighted.mean(t$event,t$ipw_trim)})
  allrisk[[paste(cohort,g)]]<-data.frame(cohort,childhood_symptom=g,first_interval_risk=h[1],second_interval_risk=h[2],cumulative_risk=1-prod(1-h),baseline_n=sum(b$childhood_symptom==g))
 }
 cat(cohort,'people',length(unique(z$id)),'events',sum(z$event),'\n')
}
write.csv(do.call(rbind,allres),file.path(out,'NCDS_BCS70_adult_interval_HR.csv'),row.names=FALSE)
write.csv(do.call(rbind,allflow),file.path(out,'NCDS_BCS70_adult_sample_flow.csv'),row.names=FALSE)
write.csv(do.call(rbind,allrisk),file.path(out,'NCDS_BCS70_adult_absolute_risk.csv'),row.names=FALSE)
write.csv(do.call(rbind,phres),file.path(out,'NCDS_BCS70_childhood_effect_time_check.csv'),row.names=FALSE)
