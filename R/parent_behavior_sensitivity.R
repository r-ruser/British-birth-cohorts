out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
library(survey);library(sandwich)
options(survey.lonely.psu='adjust')
# Load only helper definitions, without rerunning the main analysis.
ee<-parse('parent_behavior_models.R')
for(e in ee) if(is.call(e) && identical(e[[1]],as.name('<-')) && is.call(e[[3]]) && identical(e[[3]][[1]],as.name('function'))) eval(e)
ans<-list();k<-0L
for(cohort in c('MCS','NCDS','BCS70')) {
 d<-readRDS(file.path(out,paste0(cohort,'_parent_behavior_dataset.rds')))
 d<-d[complete.cases(d[c('mother_smoke','father_smoke')]),]
 for(exposure in c('mother_smoke','father_smoke','any_parent_smoke')) {
  vv<-c('male','ses','pregnancy_continued_smoke',exposure)
  if(exposure=='mother_smoke') vv<-c(vv,'father_smoke')
  if(exposure=='father_smoke') vv<-c(vv,'mother_smoke')
  level<-if(exposure=='any_parent_smoke') 'both_parent_status_known' else 'mutually_adjusted_parents'
  if(cohort=='MCS') {
   for(endpoint in c('first_asthma_report','wheeze_age11_after_age3_wheeze')) {
    if(endpoint=='first_asthma_report') {
     z<-mcperiods(d,'ever_asthma3',c(5,7,11,14),'ever_asthma',vv)
     form<-as.formula(paste('event ~ interval +',paste(vv,collapse='+'),'+ offset(log(duration))'))
     fam<-quasibinomial(link='cloglog');scale<-'HR'
    } else {
     z<-d[!is.na(d$wheeze3)&d$wheeze3==1&!is.na(d$wheeze11)&complete.cases(d[c(vv,'weight_age11','family','design_stratum')]),]
     z$event<-z$wheeze11;z$weight<-z$weight_age11;z$stratum<-z$design_stratum
     form<-reformulate(vv,'event');fam<-quasipoisson(link='log');scale<-'RR'
    }
    des<-svydesign(ids=~family,strata=~stratum,weights=~weight,data=z,nest=TRUE)
    f<-svyglm(form,design=des,family=fam);k<-k+1L
    ans[[k]]<-record(f,vcov(f),cohort,endpoint,level,exposure,length(unique(z$id)),sum(z$event),nrow(z),scale)
   }
  } else {
   ages<-if(cohort=='NCDS') c(50,55,62) else c(38,42,46);ys<-paste0('y',ages)
   a<-d[complete.cases(d[vv]),]
   response<-function(a,y) {a$response<-as.integer(!is.na(a[[y]]));f<-glm(reformulate(vv,'response'),data=a,family=binomial());pp<-pmin(pmax(predict(f,type='response'),.05),.95);mean(a$response)/pp}
   a$w0<-response(a,ys[1]);b<-a[!is.na(a[[ys[1]]])&a[[ys[1]]]==0,]
   b$w1<-b$w0*response(b,ys[2]);c<-b[!is.na(b[[ys[2]]])&b[[ys[2]]]==0,];c$w2<-c$w1*response(c,ys[3])
   r1<-b[!is.na(b[[ys[2]]]),c('id',vv,'w1',ys[2])];names(r1)[names(r1)=='w1']<-'weight';names(r1)[names(r1)==ys[2]]<-'event';r1$interval<-paste0(ages[1],'-',ages[2]);r1$duration<-diff(ages)[1]
   r2<-c[!is.na(c[[ys[3]]]),c('id',vv,'w2',ys[3])];names(r2)[names(r2)=='w2']<-'weight';names(r2)[names(r2)==ys[3]]<-'event';r2$interval<-paste0(ages[2],'-',ages[3]);r2$duration<-diff(ages)[2]
   z<-rbind(r1,r2);z$interval<-factor(z$interval);qq<-quantile(z$weight,c(.01,.99));z$weight<-pmin(pmax(z$weight,qq[1]),qq[2])
   form<-as.formula(paste('event ~ interval +',paste(vv,collapse='+'),'+ offset(log(duration))'))
   f<-glm(form,data=z,weights=weight,family=quasibinomial(link='cloglog'));V<-vcovCL(f,cluster=z$id,type='HC0');k<-k+1L
   ans[[k]]<-record(f,V,cohort,'adult_first_rereport',level,exposure,length(unique(z$id)),sum(z$event),nrow(z),'HR')
  }
 }
 cat('sensitivity ',cohort,'\n');flush.console()
}
r<-do.call(rbind,ans);r$q_BH<-NA_real_
for(f in unique(paste(r$cohort,r$endpoint))) {ix<-which(paste(r$cohort,r$endpoint)==f);r$q_BH[ix]<-p.adjust(r$p[ix],'BH')}
write.csv(r,file.path(out,'parent_behavior_complete_parents_sensitivity.csv'),row.names=FALSE)
