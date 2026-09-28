out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
fit<-readRDS(file.path(out,'MCS_first_asthma_report_weighted_model.rds'))
z<-readRDS(file.path(out,'MCS_asthma_first_report_intervals.rds'))
vv<-c('male','low_birthweight','low_income','room_smoke','respondent_asthma','atopy3','bmi3_z','ever_wheeze3')
z<-z[complete.cases(z[c('event','weight','family','stratum',vv)]),]
z$interval<-factor(z$interval,levels=c('3-5','5-7','7-11','11-14'))
cat('fit rows',length(fitted(fit)),'data rows',nrow(z),'\n')
print(head(cbind(y=z$event,fit=fitted(fit),p1=predict(fit,type='response'),p2=predict(fit,newdata=z,type='response'),duration=z$duration)))
cat('weighted obs',weighted.mean(z$event,z$weight),'fitted',weighted.mean(fitted(fit),z$weight),'predict',weighted.mean(predict(fit,newdata=z,type='response'),z$weight),'\n')
