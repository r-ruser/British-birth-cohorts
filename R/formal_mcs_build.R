root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(file.path(root,'MCS_full.rds'))
fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
getd<-function(fn){i<-which(fns==fn & x$source_inventory$format=='.dta');stopifnot(length(i)==1L);x$datasets[[i]]}
num<-function(z) suppressWarnings(as.numeric(z))
yn<-function(z) ifelse(num(z)==1,1L,ifelse(num(z)==2,0L,NA_integer_))
makeid<-function(d,p) paste(as.character(d$MCSID),num(d[[paste0(p,'CNUM00')]]),sep=':')

d<-getd('mcs1_cm_derived.dta')
a<-data.frame(id=makeid(d,'A'),family=as.character(d$MCSID),birthweight=num(d$ADBWGT00),gest_days=num(d$ADGEST00))
a$birthweight[!is.finite(a$birthweight)|a$birthweight<.5|a$birthweight>7]<-NA_real_
a$gest_days[!is.finite(a$gest_days)|a$gest_days<160|a$gest_days>330]<-NA_real_
stopifnot(!anyDuplicated(a$id))

d<-getd('mcs1_hhgrid.dta')
b<-data.frame(id=makeid(d,'A'),sex=num(d$AHCSEX00))
b<-b[b$sex %in% c(1,2),];stopifnot(!anyDuplicated(b$id))
a<-merge(a,b,by='id',all.x=TRUE,sort=FALSE)

d<-getd('mcs1_family_derived.dta')
b<-data.frame(family=as.character(d$MCSID),income_q=num(d$AOECDUK0),baseline_weight=num(d$ADWGT200),stratum=num(d$ADPTY200))
stopifnot(!anyDuplicated(b$family));a<-merge(a,b,by='family',all.x=TRUE,sort=FALSE)

d<-getd('mcs1_parent_interview.dta')
b<-data.frame(family=as.character(d$MCSID),pnum=num(d$APNUM00),room_smoke=yn(d$APSMKR00),respondent_asthma=yn(d$APASTH00))
b<-b[b$pnum==1,];stopifnot(!anyDuplicated(b$family));b$pnum<-NULL
a<-merge(a,b,by='family',all.x=TRUE,sort=FALSE)

d<-getd('mcs2_parent_cm_interview.dta')
b<-data.frame(id=makeid(d,'B'),ever_asthma3=yn(d$BPASMA00),ever_wheeze3=yn(d$BPWHEE00),atopy3=yn(d$BPECHA00))
ew<-num(d$BPWHEE00);cw<-num(d$BPWHLY00)
b$wheeze3<-ifelse(ew==2,0L,ifelse(ew==1 & cw==1,1L,ifelse(ew==1 & cw==2,0L,NA_integer_)))
b$score<-rowSums(!is.na(b[c('ever_asthma3','ever_wheeze3','wheeze3','atopy3')]))
b<-b[order(b$id,-b$score),];b<-b[!duplicated(b$id),];b$score<-NULL

e<-getd('mcs2_cm_interview.dta')
c<-data.frame(id=makeid(e,'B'),bmi3=num(e$BCBMIN00));c<-c[!duplicated(c$id),]
b<-merge(b,c,by='id',all.x=TRUE,sort=FALSE)
a<-merge(b,a,by='id',all.x=TRUE,sort=FALSE)
rm(x);gc()

# Add later waves from the independently checked outcome extraction.
l<-readRDS('E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs/MCS_analysis_long.rds')
for(age in c(5,7,11)) {
 t<-l[l$age==age,c('id','wheeze12m','ever_asthma')]
 names(t)[-1]<-paste0(c('wheeze','ever_asthma'),age)
 stopifnot(!anyDuplicated(t$id));a<-merge(a,t,by='id',all.x=TRUE,sort=FALSE)
}
# Age 14 ever-asthma outcome from the checked source table.
root2<-readRDS(file.path(root,'MCS_full.rds'))
fns2<-basename(gsub('\\\\','/',root2$source_inventory$source_entry))
e<-root2$datasets[[which(fns2=='mcs6_parent_cm_interview.dta')]]
t<-data.frame(id=makeid(e,'F'),ever_asthma14=yn(e$FPASMA00))
t$score<-!is.na(t$ever_asthma14);t<-t[order(t$id,-t$score),];t<-t[!duplicated(t$id),c('id','ever_asthma14')]
a<-merge(a,t,by='id',all.x=TRUE,sort=FALSE)
e<-root2$datasets[[which(fns2=='mcs_longitudinal_cm_file.dta')]]
t<-data.frame(id=paste(as.character(e$mcsid),num(e$cnum),sep=':'),weight2=num(e$weight2),design_stratum=num(e$pttype2),
    weight_age5=num(e$covwt2),weight_age7=num(e$dovwt2),weight_age11=num(e$eovwt2),weight_age14=num(e$fovwt2),
    age11_nonresponse_weight=num(e$enrespwt))
stopifnot(!anyDuplicated(t$id));a<-merge(a,t,by='id',all.x=TRUE,sort=FALSE)
rm(root2,l);gc()

# Preserve missing as missing, including out-of-range codes.
a$income_q[!a$income_q %in% 1:5]<-NA_real_
a$baseline_weight[!is.finite(a$baseline_weight)|a$baseline_weight<=0]<-NA_real_
a$weight2[!is.finite(a$weight2)|a$weight2<=0]<-NA_real_
a$age11_nonresponse_weight[!is.finite(a$age11_nonresponse_weight)|a$age11_nonresponse_weight<=0]<-NA_real_
for(v in c('weight_age5','weight_age7','weight_age11','weight_age14')) a[[v]][!is.finite(a[[v]])|a[[v]]<=0]<-NA_real_
a$bmi3[!is.finite(a$bmi3)|a$bmi3<8|a$bmi3>60]<-NA_real_
a$male<-ifelse(a$sex==1,1L,ifelse(a$sex==2,0L,NA_integer_))
a$low_birthweight<-ifelse(is.na(a$birthweight),NA_integer_,as.integer(a$birthweight<2.5))
a$low_income<-ifelse(is.na(a$income_q),NA_integer_,as.integer(a$income_q<=2))
bi<-quantile(a$bmi3,c(.01,.99),na.rm=TRUE)
a$bmi3[a$bmi3<bi[1]|a$bmi3>bi[2]]<-NA_real_
a$bmi3_z<-(a$bmi3-mean(a$bmi3,na.rm=TRUE))/sd(a$bmi3,na.rm=TRUE)
a$family[is.na(a$family)]<-sub(':.*$','',a$id[is.na(a$family)])

flow<-data.frame(stage=c('age3_records','ever_asthma3_no','ever_wheeze3_no','wheeze3_yes','age11_wheeze_valid','complete_predictors'),
 n=c(nrow(a),sum(a$ever_asthma3==0,na.rm=TRUE),sum(a$ever_wheeze3==0,na.rm=TRUE),sum(a$wheeze3==1,na.rm=TRUE),sum(!is.na(a$wheeze11)),sum(complete.cases(a[c('male','low_birthweight','low_income','room_smoke','respondent_asthma','atopy3','bmi3_z')]))))
write.csv(flow,file.path(out,'MCS_sample_flow.csv'),row.names=FALSE)
vars<-c('male','birthweight','low_birthweight','income_q','low_income','room_smoke','respondent_asthma','atopy3','bmi3','bmi3_z','ever_asthma3','ever_wheeze3','wheeze3','wheeze5','wheeze7','wheeze11','ever_asthma5','ever_asthma7','ever_asthma11','ever_asthma14','baseline_weight','weight2','weight_age5','weight_age7','weight_age11','weight_age14','age11_nonresponse_weight')
miss<-data.frame(variable=vars,n_total=nrow(a),n_missing=vapply(a[vars],function(z) sum(is.na(z)),0L))
miss$missing_percent<-round(100*miss$n_missing/miss$n_total,2)
write.csv(miss,file.path(out,'MCS_variable_missingness.csv'),row.names=FALSE)
saveRDS(a,file.path(out,'MCS_formal_analysis_dataset.rds'),compress='xz')
cat('Built',nrow(a),'age3 children; birth link',sum(!is.na(a$birthweight)),'income',sum(!is.na(a$income_q)),'age11 outcome',sum(!is.na(a$wheeze11)),'\n')
