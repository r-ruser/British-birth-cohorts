root<-'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out<-'E:/cohort_rds_stage_20260925/life_course_v2_20260928'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,'source_code_audit.txt'))
raw<-list()
for(cohort in c('NCDS','BCS70','MCS')) {
 cat('\nCOHORT ',cohort,'\n');flush.console()
 x<-readRDS(file.path(root,paste0(cohort,'_full.rds')))
 fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
 cat('FILES ',paste(fns,collapse='; '),'\n')
 targets<-switch(cohort,NCDS=c('ncds0123.dta','ncds4.dta','ncds5cmi.dta','ncds37yearsurveydata.dta','ncds6_v2.dta','ncds_2008_followup.dta','ncds_2013_flatfile.dta','ncds10_age62_main_interview.dta'),
  BCS70=c('sn3723.dta','bcs96x.dta','bcs2000.dta','bcs_2004_followup.dta','bcs_2008_followup.dta','bcs70_2012_flatfile.dta','bcs_age46_main.dta'),
  MCS=c('mcs_longitudinal_cm_file.dta'))
 for(fn in targets) {
  ii<-which(fns==fn);if(length(ii)!=1)next
  d<-x$datasets[[ii]];cat('\nFILE ',fn,'ROWS ',nrow(d),'\n')
  labs<-vapply(d,function(z){l<-attr(z,'label');if(is.null(l))'' else as.character(l)[1]},'')
  pat<-if(cohort=='MCS') 'weight|strat|outcome|death|dead|died' else 'asthma|wheez|long.?term|chronic|death|dead|died|n1306|b7_20|n1305|b7_1|^ncdsid|^bcsid'
  jj<-grep(pat,paste(names(d),labs),ignore.case=TRUE)
  kept<-unique(c(grep('^ncdsid$|^bcsid$|^mcsid$|^cnum$',names(d),ignore.case=TRUE),jj))
  raw[[paste(cohort,fn,sep=':')]]<-d[kept]
  for(j in jj) {
   z<-d[[j]];vl<-attr(z,'labels');tb<-sort(table(as.character(z),useNA='ifany'),decreasing=TRUE)
   cat(names(d)[j],' :: ',labs[j],' :: ',paste(names(vl),vl,sep=':',collapse='|'),' :: FREQ ',paste(names(head(tb,22)),head(tb,22),sep=':',collapse='|'),'\n')
  }
 }
 rm(x);gc()
}
sink()
saveRDS(raw,file.path(out,'selected_source_tables.rds'),compress='gzip')
cat('saved ',length(raw),'tables\n')
