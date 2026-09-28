root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/trajectory_key_and_gate_audit.txt'
sink(out)
for (cohort in c('MCS','NextSteps','BCS70','NCDS')) {
  cat('\n###',cohort,'\n')
  x <- readRDS(file.path(root,paste0(cohort,'_full.rds')))
  files <- switch(cohort,
    MCS=c('mcs2_parent_cm_interview.dta','mcs6_parent_cm_interview.dta'),
    NextSteps=c('ns8_2015_main_interview.dta','ns9_2022_main_interview.dta'),
    BCS70=c('sn3723.dta','bcs96x.dta','bcs2000.dta','bcs_2004_followup.dta'),
    NCDS=c('ncds6_v2.dta'))
  for (filename in files) {
    ix <- which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==filename &
                  x$source_inventory$format=='.dta')
    if (length(ix)!=1L) next
    d <- x$datasets[[ix]]
    cat('\nFILE',filename,'n=',nrow(d),'\n')
    vars <- names(d)
    labs <- vapply(d,function(v) if(is.null(attr(v,'label'))) '' else attr(v,'label'),character(1))
    match <- grep('cnum|person|child number|cohort member|long.?stand|respir|wheez|asthma|tuber',
                  paste(vars,labs),ignore.case=TRUE)
    if (length(match)>80) match <- match[seq_len(80)]
    for (j in match) cat(vars[j],' :: ',labs[j],'\n',sep='')
  }
  rm(x)
  gc()
}
sink()
cat(out,'\n')
