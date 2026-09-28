root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/formal_predictor_details.txt'
sink(out)
for(cohort in c('MCS','NCDS','BCS70','NextSteps')) {
 cat('\n##',cohort,'\n')
 x<-readRDS(file.path(root,paste0(cohort,'_full.rds')))
 fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
 cat('WEIGHT/RESPONSE FILES:',paste(fns[grepl('weight|response|survey_design',fns,ignore.case=TRUE)],collapse='; '),'\n')
 target<-switch(cohort,
  MCS=c('mcs1_parent_interview.dta','mcs1_parent_cm_interview.dta','mcs1_hhgrid.dta','mcs1_cm_derived.dta','mcs1_family_derived.dta','mcs2_parent_cm_interview.dta'),
  NCDS=c('ncds0123.dta','ncds_2013_derived.dta'),
  BCS70=c('sn3723.dta','bcs70_2012_derived.dta'),
  NextSteps=c('ns8_2015_derived.dta','ns8_2015_main_interview.dta'))
 for(fn in target) {
  i<-which(fns==fn);if(length(i)!=1) next
  d<-x$datasets[[i]];cat('\n',fn,'n',nrow(d),'\n')
  labs<-vapply(d,function(z) if(is.null(attr(z,'label'))) '' else as.character(attr(z,'label'))[1],'')
  pat<-if(cohort=='MCS') 'smok|asth|wheez|allerg|eczema|sex|birth weight|gestat|income quint|occupation|social class|weight|strat|cluster' else 'smok|asth|wheez|sex|birth weight|social class|income|weight|strat|cluster'
  hit<-grep(pat,paste(names(d),labs),ignore.case=TRUE)
  if(length(hit)>250) hit<-hit[1:250]
  for(j in hit) cat(names(d)[j],'::',labs[j],'\n')
 }
 rm(x);gc()
}
sink();cat(out,'\n')
