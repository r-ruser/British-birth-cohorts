root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/trajectory_variable_qc.csv'
spec <- data.frame(
  cohort=c(rep('NCDS',8),rep('BCS70',8),rep('MCS',6),rep('NextSteps',2)),
  file=c('ncds0123.dta','ncds0123.dta','ncds0123.dta','ncds4.dta',
         'ncds5cmi.dta','ncds6_v2.dta','ncds_2008_followup.dta','ncds_2013_flatfile.dta',
         'sn3723.dta','bcs7016x.dta','bcs96x.dta','bcs2000.dta',
         'bcs_2004_followup.dta','bcs_2008_followup.dta','bcs70_2012_flatfile.dta','bcs_age46_main.dta',
         'mcs2_parent_cm_interview.dta','mcs3_parent_cm_interview.dta','mcs4_parent_cm_interview.dta',
         'mcs5_parent_cm_interview.dta','mcs6_parent_cm_interview.dta','mcs8_23y_cm_survey.dta',
         'ns8_2015_main_interview.dta','ns9_2022_main_interview.dta'),
  variable=c('n259','n1305','n2617','n5770','n504021','cl112m4','N8KHPB01','N9KHPB01',
             'meb4_26','ha4_6','b960514','cl112m4','b7asth2m','b8khpb01','B9KHPB02','B10KHPB01',
             'BPASMA00','CPASMA00','DPASMA00','EPASMA00','FPASMA00','hsclsm30_m',
             'W8LOLR0Z','W9DLOILCOND0A'),
  age=c(7,11,16,23,33,42,50,55,10,16,26,30,34,38,42,46,3,5,7,11,14,23,25,32),
  stringsAsFactors=FALSE
)
rows <- list()
k <- 0L
for (cohort in unique(spec$cohort)) {
  cat('QC',cohort,'\n')
  bundle <- readRDS(file.path(root,paste0(cohort,'_full.rds')))
  inv <- bundle$source_inventory
  for (s in which(spec$cohort==cohort)) {
    ix <- which(basename(gsub('\\\\','/',inv$source_entry))==spec$file[s] & inv$format=='.dta')
    stopifnot(length(ix)==1L)
    data <- bundle$datasets[[ix]]
    v <- spec$variable[s]
    stopifnot(v %in% names(data))
    z <- data[[v]]
    counts <- sort(table(as.character(z),useNA='always'),decreasing=TRUE)
    labels <- attr(z,'labels')
    id <- grep('^(ncdsid|bcsid|mcsid|cnum|nsid)$',names(data),ignore.case=TRUE,value=TRUE)
    k <- k+1L
    rows[[k]] <- data.frame(cohort=cohort,age=spec$age[s],file=spec$file[s],
                            variable=v,n=nrow(data),id_columns=paste(id,collapse='|'),
                            variable_label=if(is.null(attr(z,'label'))) '' else attr(z,'label'),
                            value_labels=if(is.null(labels)) '' else paste(names(labels),labels,sep=':',collapse='|'),
                            values=paste(names(counts),as.integer(counts),sep=':',collapse='|'),
                            stringsAsFactors=FALSE)
  }
  rm(bundle)
  gc()
}
write.csv(do.call(rbind,rows),out,row.names=FALSE,fileEncoding='UTF-8')
cat('Wrote',out,'\n')
