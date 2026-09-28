root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/design_variables.txt'
sink(out)
for(cohort in c('MCS','NCDS','BCS70','NextSteps')) {
 x<-readRDS(file.path(root,paste0(cohort,'_full.rds'))); fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
 cat('\n##',cohort,'\n')
 target<-switch(cohort,MCS=c('mcs1_cm_derived.dta','mcs1_hhgrid.dta','mcs1_family_derived.dta','mcs2_family_derived.dta','mcs3_family_derived.dta','mcs4_family_derived.dta','mcs5_family_derived.dta','mcs6_family_derived.dta'),NCDS=c('ncds_response.dta'),BCS70=c('sn3723.dta','bcs_age46_main.dta'),NextSteps=c('ns8_2015_derived.dta','ns9_2022_derived_variables.dta'))
 for(fn in target) {
  i<-which(fns==fn);if(length(i)!=1) next
  d<-x$datasets[[i]];cat('\n',fn,'n=',nrow(d),'\n')
  labs<-vapply(d,function(z) if(is.null(attr(z,'label'))) '' else as.character(attr(z,'label'))[1],'')
  hit<-grep('mcsid|cnum|pnum|sex|weight.*analysis|survey weight|longitudinal weight|non.?response weight|^.*WGT[0-9]|strat|cluster|^NSID$|^NCDSID$',paste(names(d),labs),ignore.case=TRUE)
  if(length(hit)>100) hit<-hit[1:100]
  for(j in hit) cat(names(d)[j],'::',labs[j],'\n')
 }
 rm(x);gc()
}
sink();cat(out,'\n')
