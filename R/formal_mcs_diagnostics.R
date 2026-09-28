out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
options(survey.lonely.psu='adjust')
library(survey)
z<-readRDS(file.path(out,'MCS_asthma_first_report_intervals.rds'))
vv<-c('male','low_birthweight','low_income','room_smoke','respondent_asthma','atopy3','bmi3_z','ever_wheeze3')
z<-z[complete.cases(z[c('event','weight','family','stratum',vv)]),]
z$interval<-factor(z$interval,levels=c('3-5','5-7','7-11','11-14'))
des<-svydesign(ids=~family,strata=~stratum,weights=~weight,data=z,nest=TRUE)
base<-svyglm(event~interval+male+low_birthweight+low_income+room_smoke+respondent_asthma+atopy3+bmi3_z+ever_wheeze3+offset(log(duration)),design=des,family=quasibinomial(link='cloglog'))
rows<-list()
for(v in c('male','atopy3','ever_wheeze3','bmi3_z')) {
 f<-update(formula(base),as.formula(paste('. ~ . + interval:',v)))
 alt<-svyglm(f,design=des,family=quasibinomial(link='cloglog'))
 term<-paste0('interval:',v)
 test<-regTermTest(alt,term)
 rows[[v]]<-data.frame(predictor=v,term=paste(deparse(test$call),collapse=''),Wald_F=as.numeric(test$Ftest),df=as.numeric(test$df),p=as.numeric(test$p))
}
write.csv(do.call(rbind,rows),file.path(out,'MCS_HR_time_varying_effect_checks.csv'),row.names=FALSE)

alt<-svyglm(update(formula(base),. ~ . + interval:ever_wheeze3),design=des,family=quasibinomial(link='cloglog'))
b<-coef(alt);V<-vcov(alt)
tb<-list()
for(iv in levels(z$interval)) {
  L<-setNames(rep(0,length(b)),names(b));L['ever_wheeze3']<-1
  if(iv!=levels(z$interval)[1]) {nm<-paste0('interval',iv,':ever_wheeze3');if(!nm%in%names(L)) nm<-paste0('ever_wheeze3:interval',iv);L[nm]<-1}
  est<-sum(L*b);se<-sqrt(as.numeric(t(L)%*%V%*%L))
  tb[[iv]]<-data.frame(interval=iv,HR=exp(est),CI_low=exp(est-1.96*se),CI_high=exp(est+1.96*se))
}
write.csv(do.call(rbind,tb),file.path(out,'MCS_baseline_wheeze_interval_specific_HR.csv'),row.names=FALSE)

# Weighted observed vs predicted interval event probability.
z$pred<-as.numeric(fitted(base))
diag<-do.call(rbind,lapply(levels(z$interval),function(iv){a<-z[z$interval==iv,];data.frame(interval=iv,n=nrow(a),events=sum(a$event),weighted_observed=weighted.mean(a$event,a$weight),weighted_predicted=weighted.mean(a$pred,a$weight))}))
write.csv(diag,file.path(out,'MCS_HR_interval_fit_check.csv'),row.names=FALSE)
