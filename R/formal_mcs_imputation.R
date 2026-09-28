out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
library(mice);library(survey)
options(survey.lonely.psu='adjust')
set.seed(20260928)
pool_terms<-function(fits,terms,model,scale) {
 m<-length(fits);B<-sapply(fits,function(f) coef(f)[terms]);U<-sapply(fits,function(f) diag(vcov(f))[terms])
 if(is.vector(B)) B<-matrix(B,nrow=length(terms));if(is.vector(U)) U<-matrix(U,nrow=length(terms))
 q<-rowMeans(B);ubar<-rowMeans(U);between<-apply(B,1,var);total<-ubar+(1+1/m)*between
 df<-ifelse(between==0,Inf,(m-1)*(1+ubar/((1+1/m)*between))^2)
 crit<-qt(.975,df);tval<-q/sqrt(total)
 data.frame(model,predictor=terms,estimate=q,ratio=exp(q),CI_low=exp(q-crit*sqrt(total)),CI_high=exp(q+crit*sqrt(total)),p=2*pt(-abs(tval),df),imputations=m,scale)
}
run_mice<-function(wide,impute_vars,aux_vars) {
 tab<-wide[c(impute_vars,aux_vars)]
 for(v in impute_vars) if(v!='bmi3_z') tab[[v]]<-as.factor(tab[[v]])
 if('design_stratum'%in%names(tab)) tab$design_stratum<-as.factor(tab$design_stratum)
 meth<-rep('',ncol(tab));names(meth)<-names(tab)
 meth[impute_vars]<-ifelse(impute_vars=='bmi3_z','pmm','logreg')
 pm<-make.predictorMatrix(tab);diag(pm)<-0
 for(v in aux_vars) pm[v,]<-0
 # Avoid imputation of identification/design auxiliaries and allow them to inform missing predictors.
 im<-mice(tab,m=5,maxit=5,method=meth,predictorMatrix=pm,printFlag=FALSE,seed=20260928)
 lapply(seq_len(5),function(i){a<-complete(im,i);for(v in impute_vars) if(v!='bmi3_z') a[[v]]<-as.numeric(as.character(a[[v]]));a})
}

covs<-c('male','low_birthweight','low_income','room_smoke','respondent_asthma','atopy3','bmi3_z','ever_wheeze3')
z<-readRDS(file.path(out,'MCS_asthma_first_report_intervals.rds'))
wide<-z[!duplicated(z$id),c('id','family','stratum',covs)]
event_summary<-aggregate(event~id,z,max);names(event_summary)[2]<-'ever_event'
obs<-aggregate(event~id,z,length);names(obs)[2]<-'observed_intervals'
wide<-merge(wide,event_summary,by='id',sort=FALSE);wide<-merge(wide,obs,by='id',sort=FALSE)
wide$design_stratum<-wide$stratum;wide$stratum<-NULL
imp<-run_mice(wide,covs,c('ever_event','observed_intervals','design_stratum'))
fits<-list()
for(i in seq_along(imp)) {
 a<-merge(z[c('id','family','stratum','interval','duration','event','weight')],cbind(id=wide$id,imp[[i]][covs]),by='id',sort=FALSE)
 a$interval<-factor(a$interval,levels=c('3-5','5-7','7-11','11-14'))
 des<-svydesign(ids=~family,strata=~stratum,weights=~weight,data=a,nest=TRUE)
 fits[[i]]<-svyglm(event~interval+male+low_birthweight+low_income+room_smoke+respondent_asthma+atopy3+bmi3_z+ever_wheeze3+offset(log(duration)),design=des,family=quasibinomial(link='cloglog'))
}
r1<-pool_terms(fits,covs,'MCS_first_asthma_report','HR')
cat('asthma MI',nrow(wide),'people',sum(z$event),'events\n')
write.csv(r1,file.path(out,'MCS_asthma_MI_partial.csv'),row.names=FALSE)

d<-readRDS(file.path(out,'MCS_formal_analysis_dataset.rds'))
pro<-d[!is.na(d$wheeze3) & d$wheeze3==1 & !is.na(d$wheeze11) & !is.na(d$weight_age11) & !is.na(d$family) & !is.na(d$design_stratum),]
pv<-c(setdiff(covs,'ever_wheeze3'),'ever_asthma3')
pro$design_stratum<-pro$design_stratum
imp2<-run_mice(pro,pv,c('wheeze11','design_stratum'))
fits2<-list()
for(i in seq_along(imp2)) {
 a<-pro[c('id','family','design_stratum','weight_age11','wheeze11')]
 a[pv]<-imp2[[i]][pv]
 des<-svydesign(ids=~family,strata=~design_stratum,weights=~weight_age11,data=a,nest=TRUE)
 fits2[[i]]<-svyglm(wheeze11~male+low_birthweight+low_income+room_smoke+respondent_asthma+atopy3+bmi3_z+ever_asthma3,design=des,family=quasipoisson(link='log'))
}
r2<-pool_terms(fits2,pv,'MCS_age3_wheeze_to_age11','RR')
cat('prognosis MI',nrow(pro),'people',sum(pro$wheeze11),'events\n')
write.csv(rbind(r1,r2),file.path(out,'MCS_multiple_imputation_sensitivity.csv'),row.names=FALSE)
