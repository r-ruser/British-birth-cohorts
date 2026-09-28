root<-'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
x<-readRDS(file.path(root,'NextSteps_full.rds'))
fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
getd<-function(fn){i<-which(fns==fn);stopifnot(length(i)==1);x$datasets[[i]]}
num<-function(z) suppressWarnings(as.numeric(z))
yn<-function(z) ifelse(num(z)==1,1L,ifelse(num(z)==2,0L,NA_integer_))
a<-getd('ns8_2015_main_interview.dta');n25<-nrow(a)
d<-data.frame(id=as.character(a$NSID),sex=num(a$W8CMSEX),illness25=yn(a$W8LOIL),impact25=yn(a$W8LOLR0Z))
a<-getd('ns8_2015_derived.dta');b<-data.frame(id=as.character(a$NSID),bmi25=num(a$W8DBMI));d<-merge(d,b,by='id',all.x=TRUE,sort=FALSE)
a<-getd('ns9_2022_main_interview.dta');n32<-nrow(a);b<-data.frame(id=as.character(a$NSID),impact32=yn(a$W9LOLRAD),resp_cond32=yn(a$W9DLOILCOND0A))
d<-merge(d,b,by='id',all=TRUE,sort=FALSE)
d$bmi25[!is.finite(d$bmi25)|d$bmi25<15|d$bmi25>60]<-NA_real_
d$bmi_cat<-cut(d$bmi25,c(0,25,30,Inf),right=FALSE,labels=c('under_25','25_to_29.9','30_plus'))
flow<-data.frame(stage=c('age25_records','age32_records','age32_respiratory_condition_question_valid','age32_valid_with_age25_BMI','age32_valid_with_age25_sex'),
 n=c(n25,n32,sum(!is.na(d$resp_cond32)),sum(!is.na(d$resp_cond32)&!is.na(d$bmi_cat)),sum(!is.na(d$resp_cond32)&d$sex%in%c(1,2))))
write.csv(flow,file.path(out,'NextSteps_sample_flow.csv'),row.names=FALSE)
rows<-list()
for(var in c('bmi_cat','sex','illness25')) {
 z<-d[!is.na(d$resp_cond32)&!is.na(d[[var]]),]
 for(level in unique(z[[var]])) {
  q<-z[z[[var]]==level,];rows[[length(rows)+1L]]<-data.frame(group_variable=var,group=as.character(level),n=nrow(q),respiratory_yes=sum(q$resp_cond32),proportion=mean(q$resp_cond32))
 }
}
write.csv(do.call(rbind,rows),file.path(out,'NextSteps_age32_conditional_respiratory_descriptive.csv'),row.names=FALSE)
cat('age32 respiratory condition',sum(d$resp_cond32==1,na.rm=TRUE),'/',sum(!is.na(d$resp_cond32)),'\n')
