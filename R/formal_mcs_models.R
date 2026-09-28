out <- 'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
d <- readRDS(file.path(out,'MCS_formal_analysis_dataset.rds'))
options(survey.lonely.psu='adjust')
library(survey)
library(pROC)

basevars <- c('male','low_birthweight','low_income','room_smoke','respondent_asthma','atopy3','bmi3_z')
make_periods <- function(data,basevar,ages,stem) {
  z <- data[data[[basevar]]==0 & !is.na(data[[basevar]]),]
  rows <- vector('list',nrow(z));n <- 0L
  for(i in seq_len(nrow(z))) {
    a <- z[i,,drop=FALSE];one <- list();j <- 0L;prev <- 3
    for(age in ages) {
      y <- a[[paste0(stem,age)]]
      if(is.na(y)) break
      wt <- a[[paste0('weight_age',age)]]
      if(is.na(wt)) break
      j<-j+1L
      one[[j]]<-data.frame(id=a$id,family=a$family,stratum=a$design_stratum,
        interval=paste0(prev,'-',age),duration=age-prev,event=y,weight=wt,
        a[setdiff(c(basevars,'ever_wheeze3'),basevar)],check.names=FALSE)
      if(y==1) break
      prev<-age
    }
    if(j) {n<-n+1L;rows[[n]]<-do.call(rbind,one)}
  }
  do.call(rbind,rows[seq_len(n)])
}
hrfit <- function(z,variables,label) {
  z <- z[complete.cases(z[c('event','weight','family','stratum',variables)]),]
  z$interval<-factor(z$interval,levels=unique(z$interval))
  stopifnot(all(z$weight>0),all(z$duration>0))
  des <- svydesign(ids=~family,strata=~stratum,weights=~weight,data=z,nest=TRUE)
  form <- as.formula(paste('event ~ interval +',paste(variables,collapse='+'),'+ offset(log(duration))'))
  fit <- svyglm(form,design=des,family=quasibinomial(link='cloglog'))
  b <- coef(fit);se <- sqrt(diag(vcov(fit)))
  hit <- names(b)%in%variables
  result<-data.frame(model=label,predictor=names(b)[hit],HR=exp(b[hit]),CI_low=exp(b[hit]-1.96*se[hit]),CI_high=exp(b[hit]+1.96*se[hit]),p=2*pnorm(-abs(b[hit]/se[hit])),
    people=length(unique(z$id)),interval_rows=nrow(z),events=sum(z$event),stringsAsFactors=FALSE)
  flow<-aggregate(event~interval,z,function(v)c(at_risk=length(v),events=sum(v)))
  flow<-data.frame(model=label,interval=flow$interval,at_risk=flow$event[,1],events=flow$event[,2])
  saveRDS(fit,file.path(out,paste0(label,'_weighted_model.rds')),compress='xz')
  write.csv(flow,file.path(out,paste0(label,'_risk_intervals.csv')),row.names=FALSE)
  cat(label,'people',result$people[1],'events',result$events[1],'\n')
  result
}

# Main endpoint: first subsequent parent report of ever-asthma among those with no such report at age 3.
asth <- make_periods(d,'ever_asthma3',c(5,7,11,14),'ever_asthma')
saveRDS(asth,file.path(out,'MCS_asthma_first_report_intervals.rds'),compress='xz')
r1<-hrfit(asth,c(basevars,'ever_wheeze3'),'MCS_first_asthma_report')

# Symptom sensitivity: first observed current wheeze after age 3 among children never reported to wheeze by age 3.
wh <- make_periods(d,'ever_wheeze3',c(5,7,11),'wheeze')
saveRDS(wh,file.path(out,'MCS_wheeze_first_report_intervals.rds'),compress='xz')
r2<-hrfit(wh,basevars,'MCS_first_wheeze_report')
write.csv(rbind(r1,r2),file.path(out,'MCS_weighted_interval_HR.csv'),row.names=FALSE)

# Adverse prognosis: age-11 wheeze among children with current wheeze at age 3.
pro<-d[d$wheeze3==1 & !is.na(d$wheeze11),]
pro$weight<-pro$weight_age11
pv<-c(basevars,'ever_asthma3')
pro<-pro[complete.cases(pro[c('wheeze11','weight','family','design_stratum',pv)]),]
des<-svydesign(ids=~family,strata=~design_stratum,weights=~weight,data=pro,nest=TRUE)
form<-as.formula(paste('wheeze11 ~',paste(pv,collapse='+')))
fit<-svyglm(form,design=des,family=quasipoisson(link='log'))
b<-coef(fit);se<-sqrt(diag(vcov(fit)));j<-names(b)%in%pv
rr<-data.frame(predictor=names(b)[j],RR=exp(b[j]),CI_low=exp(b[j]-1.96*se[j]),CI_high=exp(b[j]+1.96*se[j]),p=2*pnorm(-abs(b[j]/se[j])),people=nrow(pro),events=sum(pro$wheeze11))
write.csv(rr,file.path(out,'MCS_age3_wheeze_to_age11_weighted_RR.csv'),row.names=FALSE)
saveRDS(fit,file.path(out,'MCS_age3_wheeze_to_age11_weighted_model.rds'),compress='xz')
write.csv(data.frame(stage=c('age3_current_wheeze','age11_outcome_available','model_complete'),n=c(sum(d$wheeze3==1,na.rm=TRUE),sum(d$wheeze3==1 & !is.na(d$wheeze11),na.rm=TRUE),nrow(pro)),events=c(NA,sum(d$wheeze3==1 & d$wheeze11==1,na.rm=TRUE),sum(pro$wheeze11))),file.path(out,'MCS_prognosis_flow.csv'),row.names=FALSE)
cat('prognosis people',nrow(pro),'events',sum(pro$wheeze11),'\n')

# Internal 5-fold validation of the prognostic risk score, split by family.
set.seed(20260928)
fam<-unique(pro$family);fold<-sample(rep(1:5,length.out=length(fam)))
pro$fold<-fold[match(pro$family,fam)]
pred<-rep(NA_real_,nrow(pro))
for(k in 1:5) {
  tr<-pro[pro$fold!=k,];te<-pro[pro$fold==k,]
  f<-glm(form,data=tr,family=quasibinomial(),weights=weight)
  pred[pro$fold==k]<-predict(f,newdata=te,type='response')
}
pred<-pmin(pmax(pred,1e-5),1-1e-5)
auc<-as.numeric(pROC::auc(pro$wheeze11,pred,quiet=TRUE))
brier<-mean((pro$wheeze11-pred)^2)
cal_int<-coef(glm(pro$wheeze11~offset(qlogis(pred)),family=binomial()))[1]
cal_slope<-coef(glm(pro$wheeze11~qlogis(pred),family=binomial()))[2]
metrics<-data.frame(n=nrow(pro),events=sum(pro$wheeze11),folds=5,auc=auc,brier=brier,calibration_intercept=cal_int,calibration_slope=cal_slope,
                    observed_risk=mean(pro$wheeze11),mean_predicted=mean(pred))
write.csv(metrics,file.path(out,'MCS_prognosis_internal_validation.csv'),row.names=FALSE)
set.seed(20260928)
boot<-matrix(NA_real_,300,3,dimnames=list(NULL,c('auc','brier','calibration_slope')))
ff<-unique(pro$family)
for(i in seq_len(nrow(boot))) {
 fam_i<-sample(ff,length(ff),replace=TRUE)
 ix<-unlist(lapply(fam_i,function(f) which(pro$family==f)),use.names=FALSE)
 yy<-pro$wheeze11[ix];pp<-pred[ix]
 if(length(unique(yy))<2) next
 boot[i,'auc']<-as.numeric(pROC::auc(yy,pp,quiet=TRUE))
 boot[i,'brier']<-mean((yy-pp)^2)
 boot[i,'calibration_slope']<-coef(glm(yy~qlogis(pp),family=binomial()))[2]
}
ci<-data.frame(metric=colnames(boot),estimate=c(auc,brier,cal_slope),CI_low=apply(boot,2,quantile,.025,na.rm=TRUE),CI_high=apply(boot,2,quantile,.975,na.rm=TRUE),bootstrap_replicates=300)
write.csv(ci,file.path(out,'MCS_prognosis_internal_validation_CI.csv'),row.names=FALSE)
pro$predicted_oof<-pred
pro$decile<-cut(pred,unique(quantile(pred,seq(0,1,.1),na.rm=TRUE)),include.lowest=TRUE)
cal<-aggregate(cbind(wheeze11,predicted_oof)~decile,pro,mean)
cal$n<-as.integer(table(pro$decile)[as.character(cal$decile)])
write.csv(cal,file.path(out,'MCS_prognosis_calibration.csv'),row.names=FALSE)
cat('CV AUC',auc,'Brier',brier,'cal slope',cal_slope,'\n')
