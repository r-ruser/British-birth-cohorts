root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
sink(file.path(out,'selected_source_details.txt'))
for(cohort in c('MCS','BCS70','NextSteps')) {
 cat('\n## ',cohort,'\n');x<-readRDS(file.path(root,paste0(cohort,'_full.rds')))
 fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
 cat('FILES ',paste(unique(fns),collapse='; '),'\n')
 for(i in seq_along(x$datasets)) {
  d<-x$datasets[[i]];if(!is.data.frame(d)) next
  fn<-fns[i]
  if(cohort=='MCS' && !fn%in%c('mcs1_parent_interview.dta','mcs1_parent_cm_interview.dta','mcs1_hhgrid.dta')) next
  if(cohort=='NextSteps' && !grepl('wave1|wave2|lsype.*1|lsype.*2|lsype1|lsype2',fn,ignore.case=TRUE)) next
  labs<-vapply(d,function(z){l<-attr(z,'label');if(is.null(l)) '' else as.character(l)[1]},'')
  pat<-switch(cohort,MCS='preg|sex|gender|relation|natural|biolog|APNUM|ACNUM|APSMTY|APSMEV|smok',BCS70='breast.?fed|breast.?feed|feeding|breast milk',NextSteps='parent.*smok|mother.*smok|father.*smok')
  j<-grep(pat,paste(names(d),labs),ignore.case=TRUE,perl=TRUE)
  if(length(j)) {cat('\nFILE ',fn,'n ',nrow(d),'\n')
   for(k in j) {
    z<-d[[k]];vl<-attr(z,'labels');cat(names(d)[k],'::',labs[k],'::',paste(names(vl),vl,sep=':',collapse='|'),'\n')
   }
  }
 }
 rm(x);gc()
}
sink()
