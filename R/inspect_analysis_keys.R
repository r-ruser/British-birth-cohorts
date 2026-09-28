root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/analysis_keys.txt'
spec <- list(
 NCDS=c('ncds0123.dta','ncds_2008_followup.dta','ncds_2013_flatfile.dta','ncds_2013_derived.dta','ncds10_age62_main_interview.dta'),
 BCS70=c('sn3723.dta','bcs7016x.dta','bcs_2008_followup.dta','bcs70_2012_flatfile.dta','bcs70_2012_derived.dta','bcs_age46_main.dta'),
 MCS=c('mcs2_parent_cm_interview.dta','mcs2_cm_interview.dta','mcs3_parent_cm_interview.dta','mcs3_cm_interview.dta','mcs4_parent_cm_interview.dta','mcs4_cm_interview.dta','mcs5_parent_cm_interview.dta','mcs5_cm_interview.dta','mcs6_parent_cm_interview.dta','mcs6_cm_derived.dta','mcs_core.dta'),
 NextSteps=c('ns8_2015_main_interview.dta','ns8_2015_derived.dta','ns9_2022_main_interview.dta','ns9_2022_derived_variables.dta'))
sink(out)
for (cohort in names(spec)) {
 x <- readRDS(file.path(root,paste0(cohort,'_full.rds')))
 cat('\n###',cohort,'\n')
 for (fn in spec[[cohort]]) {
  ix <- which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn & x$source_inventory$format=='.dta')
  if(length(ix)!=1L) {cat('FILE NOT FOUND',fn,'\n');next}
  d <- x$datasets[[ix]]
  cat('\n',fn,'n=',nrow(d),'\n')
  hits <- grep('^(MCSID|NSID|BCSID|NCDSID|.*CNUM00|.*PNUM00|.*BMI.*|.*ASMA.*|.*WHEE.*|.*WHLY.*|.*AGE.*)$',names(d),ignore.case=TRUE,value=TRUE)
  if(fn=='mcs_core.dta') hits <- unique(c(hits,grep('asthma|wheez|sweep|age|cohort member',vapply(d,function(z) if(is.null(attr(z,'label'))) '' else attr(z,'label'),'') ,ignore.case=TRUE,value=TRUE)))
  if (length(hits)>150) hits <- hits[1:150]
  for (v in hits) {
   z<-d[[v]]; cat(v,'::',attr(z,'label'),':: nvalid=',sum(!is.na(z)&suppressWarnings(as.numeric(z))>=0),'\n')
  }
 }
 rm(x);gc()
}
sink();cat(out,'\n')
