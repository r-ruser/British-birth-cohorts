root<-'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
num<-function(z) suppressWarnings(as.numeric(z))
yn<-function(z,yes=1,no=0) ifelse(num(z)==yes,1L,ifelse(num(z)==no,0L,NA_integer_))
getd<-function(x,fn){i<-which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn & x$source_inventory$format=='.dta');stopifnot(length(i)==1L);x$datasets[[i]]}
addout<-function(a,d,idcol,col,age,yes=1,no=0){b<-data.frame(id=as.character(d[[idcol]]),y=yn(d[[col]],yes,no));names(b)[2]<-paste0('y',age);stopifnot(!anyDuplicated(b$id));merge(a,b,by='id',all.x=TRUE,sort=FALSE)}

cat('NCDS\n');x<-readRDS(file.path(root,'NCDS_full.rds'))
d<-getd(x,'ncds0123.dta')
n259<-num(d$n259);n260<-num(d$n260);bw<-num(d$n646);sm<-num(d$n639);sc<-num(d$n492)
a<-data.frame(id=as.character(d$ncdsid),male=ifelse(num(d$n622)==1,1L,ifelse(num(d$n622)==2,0L,NA_integer_)),
 maternal_smoke=ifelse(sm %in% 2:4,1L,ifelse(sm==1,0L,NA_integer_)),
 low_social=ifelse(sc %in% 4:6,1L,ifelse(sc %in% 2:3,0L,NA_integer_)),
 low_birthweight=ifelse(bw==509,1L,ifelse(bw==508,0L,ifelse(bw>=40 & bw<=200,as.integer(bw*28.3495<2500),NA_integer_))),
 childhood_symptom=ifelse(n259==2|n260==2,1L,ifelse(n259==3 & n260==3,0L,NA_integer_)))
stopifnot(!anyDuplicated(a$id))
d<-getd(x,'ncds_2008_followup.dta');a<-addout(a,d,'NCDSID','N8KHPB01',50)
d<-getd(x,'ncds_2013_flatfile.dta');a<-addout(a,d,'NCDSID','N9KHPB01',55,1,2)
d<-getd(x,'ncds10_age62_main_interview.dta');a<-addout(a,d,'ncdsid','n10khlprb01',62)
saveRDS(a,file.path(out,'NCDS_adult_analysis_dataset.rds'),compress='xz')
rm(x,d,a);gc()

cat('BCS70\n');x<-readRDS(file.path(root,'BCS70_full.rds'))
d<-getd(x,'bcs7072a.dta');sm<-num(d$a0043b);bw<-num(d$a0278)
a<-data.frame(id=as.character(d$bcsid),male=ifelse(num(d$a0255)==1,1L,ifelse(num(d$a0255)==2,0L,NA_integer_)),
 maternal_smoke=ifelse(sm %in% 4:6,1L,ifelse(sm %in% 1:3,0L,NA_integer_)),
 low_birthweight=ifelse(bw>=500 & bw<=7000,as.integer(bw<2500),NA_integer_))
stopifnot(!anyDuplicated(a$id))
d<-getd(x,'bcs1derived.dta');sc<-num(d$BD1PSOC)
b<-data.frame(id=as.character(d$BCSID),low_social=ifelse(sc %in% 3:5,1L,ifelse(sc %in% 6:8,0L,NA_integer_)))
a<-merge(a,b,by='id',all.x=TRUE,sort=FALSE)
d<-getd(x,'sn3723.dta');b<-data.frame(id=as.character(d$bcsid),childhood_symptom=yn(d$b7_1,1,2))
a<-merge(a,b,by='id',all.x=TRUE,sort=FALSE)
d<-getd(x,'bcs_2008_followup.dta');a<-addout(a,d,'bcsid','b8khpb01',38)
d<-getd(x,'bcs70_2012_flatfile.dta');a<-addout(a,d,'BCSID','B9KHPB02',42)
d<-getd(x,'bcs_age46_main.dta');a<-addout(a,d,'BCSID','B10KHPB01',46)
d<-getd(x,'bcs_2004_followup.dta');b<-data.frame(id=as.character(d$bcsid),bmi34=num(d$bd7bmi));b$bmi34[!is.finite(b$bmi34)|b$bmi34<15|b$bmi34>60]<-NA_real_
a<-merge(a,b,by='id',all.x=TRUE,sort=FALSE)
saveRDS(a,file.path(out,'BCS70_adult_analysis_dataset.rds'),compress='xz')
rm(x,d,a);gc()
cat('done\n')
