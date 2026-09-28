root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
old <- 'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
num<-function(z) suppressWarnings(as.numeric(z))
yn<-function(z) ifelse(num(z)==1,1L,ifelse(num(z)==2,0L,NA_integer_))
anyparent<-function(m,f) ifelse(m==1|f==1,1L,ifelse(m==0&f==0,0L,NA_integer_))
getd<-function(x,fn) {i<-which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn & x$source_inventory$format=='.dta');stopifnot(length(i)==1);x$datasets[[i]]}
miss<-list();freq<-list()
for(cohort in c('MCS','NCDS','BCS70')) {
 cat('build ',cohort,'\n');flush.console()
 x<-readRDS(file.path(root,paste0(cohort,'_full.rds')))
 if(cohort=='MCS') {
  h<-getd(x,'mcs1_hhgrid.dta')
  h<-data.frame(family=as.character(h$MCSID),pnum=num(h$APNUM00),sex=num(h$AHPSEX00),relation=num(h$AHCREL00))
  h<-unique(h[h$relation%in%7:10 & h$sex%in%1:2,c('family','pnum','sex','relation')])
  stopifnot(!anyDuplicated(h[c('family','pnum')]))
  p<-getd(x,'mcs1_parent_interview.dta')
  pp<-data.frame(family=as.character(p$MCSID),pnum=num(p$APNUM00),
     current_smoke=ifelse(num(p$APSMUS0A)%in%c(2:5,95),1L,ifelse(num(p$APSMUS0A)==1,0L,NA_integer_)),
     pre=num(p$APCIPR00),change=num(p$APSMCH00),after=num(p$APCICH00),
     recent=num(p$APSMTY00),ever=num(p$APSMEV00),room=yn(p$APSMKR00))
  stopifnot(!anyDuplicated(pp[c('family','pnum')]))
  pp$pregnancy_continued_smoke<-NA_integer_
  # Derive continuation after the reported smoking change, not any prenatal exposure.
  pp$pregnancy_continued_smoke[pp$pre>=0 & pp$pre<=90 & pp$change==2 & !is.na(pp$change)]<-as.integer(pp$pre[pp$pre>=0 & pp$pre<=90 & pp$change==2 & !is.na(pp$change)]>0)
  pp$pregnancy_continued_smoke[pp$change==1 & !is.na(pp$change) & (pp$after%in%1:90 | pp$after==96)]<-1L
  pp$pregnancy_continued_smoke[pp$change==1 & !is.na(pp$change) & pp$after==0 & !is.na(pp$after)]<-0L
  pp$pregnancy_continued_smoke[pp$current_smoke==0 & !is.na(pp$current_smoke) & (pp$recent==2 | pp$ever==2) & (pp$pre<0 | is.na(pp$pre))]<-0L
  hp<-merge(h,pp,by=c('family','pnum'),all.x=TRUE,sort=FALSE)
  a<-readRDS(file.path(old,'MCS_formal_analysis_dataset.rds'))
  families<-unique(a$family)
  b<-data.frame(family=families,mother_smoke=NA_integer_,father_smoke=NA_integer_,pregnancy_continued_smoke=NA_integer_)
  for(i in seq_along(families)) {
   z<-hp[hp$family==families[i],]
   m<-z[z$sex==2,];f<-z[z$sex==1,]
   if(nrow(m)==1) b$mother_smoke[i]<-m$current_smoke
   if(nrow(f)==1) b$father_smoke[i]<-f$current_smoke
   bm<-m[m$relation==7,]
   if(nrow(bm)==1) b$pregnancy_continued_smoke[i]<-bm$pregnancy_continued_smoke
  }
  a<-merge(a,b,by='family',all.x=TRUE,sort=FALSE)
  d<-getd(x,'mcs1_parent_cm_interview.dta')
  bf<-data.frame(id=paste(as.character(d$MCSID),num(d$ACNUM00),sep=':'),pnum=num(d$APNUM00),tried=num(d$ACBFEV00),last=num(d$ACBFEA00))
  bf<-bf[bf$pnum==1,];stopifnot(!anyDuplicated(bf$id))
  bf$breastfed<-ifelse(bf$last%in%2:6,1L,ifelse(bf$last==1 | bf$tried==2,0L,NA_integer_))
  a<-merge(a,bf[c('id','breastfed')],by='id',all.x=TRUE,sort=FALSE)
  a$any_parent_smoke<-anyparent(a$mother_smoke,a$father_smoke)
  a$ses<-a$low_income
  a$behavior_age<-0.75
 } else {
  a<-readRDS(file.path(old,paste0(cohort,'_adult_analysis_dataset.rds')))
  if(cohort=='NCDS') {
   d<-getd(x,'ncds0123.dta')
   m<-num(d$n2400);f<-num(d$n2401);bf<-num(d$n222)
   b<-data.frame(id=as.character(d$ncdsid),mother_smoke=ifelse(m%in%2:8,1L,ifelse(m==1,0L,NA_integer_)),
      father_smoke=ifelse(f%in%2:8,1L,ifelse(f==1,0L,NA_integer_)),
      breastfed=ifelse(bf%in%3:4,1L,ifelse(bf==2,0L,NA_integer_)))
   age<-16
  } else {
   d<-getd(x,'sn3723.dta');m<-num(d$e9_1);f<-num(d$e11_1)
   b<-data.frame(id=as.character(d$bcsid),mother_smoke=ifelse(m%in%1:2,1L,ifelse(m==3,0L,NA_integer_)),
      father_smoke=ifelse(f%in%1:2,1L,ifelse(f==3,0L,NA_integer_)))
   b$breastfed<-NA_integer_;age<-10
  }
  stopifnot(!anyDuplicated(b$id))
  a<-merge(a,b,by='id',all.x=TRUE,sort=FALSE)
  a$any_parent_smoke<-anyparent(a$mother_smoke,a$father_smoke)
  a$pregnancy_continued_smoke<-a$maternal_smoke
  a$ses<-a$low_social;a$behavior_age<-age
 }
 stopifnot(!anyDuplicated(a$id))
 vs<-c('pregnancy_continued_smoke','mother_smoke','father_smoke','any_parent_smoke','breastfed')
 miss[[cohort]]<-data.frame(cohort,variable=vs,n=nrow(a),missing=sapply(a[vs],function(z)sum(is.na(z))))
 for(v in vs) freq[[paste(cohort,v)]]<-data.frame(cohort,variable=v,unexposed=sum(a[[v]]==0,na.rm=TRUE),exposed=sum(a[[v]]==1,na.rm=TRUE),missing=sum(is.na(a[[v]])))
 saveRDS(a,file.path(out,paste0(cohort,'_parent_behavior_dataset.rds')),compress='xz')
 rm(x);gc()
}
write.csv(do.call(rbind,miss),file.path(out,'parent_behavior_missingness.csv'),row.names=FALSE)
write.csv(do.call(rbind,freq),file.path(out,'parent_behavior_frequencies.csv'),row.names=FALSE)
