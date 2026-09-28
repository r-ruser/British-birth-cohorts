root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/selected_trajectory_values.txt'
target <- list(
  NCDS=list('ncds0123.dta'=c('n259','n260','n1305','n2617'),
            'ncds_2008_followup.dta'=c('N8KHPB01'),
            'ncds_2013_flatfile.dta'=c('N9KHPB01'),
            'ncds10_age62_main_interview.dta'=c('n10khlprb01')),
  BCS70=list('sn3723.dta'=c('b7_1','b7_2','b7_20'),
             'bcs7016x.dta'=c('ob9_1','rc3_6','ha4_5','ha4_6'),
             'bcs_2008_followup.dta'=c('b8khpb01'),
             'bcs70_2012_flatfile.dta'=c('B9KHPB02'),
             'bcs_age46_main.dta'=c('B10KHPB01')),
  MCS=list('mcs2_parent_cm_interview.dta'=c('BPWHLY00','BPASMA00','BCNUM00'),
           'mcs3_parent_cm_interview.dta'=c('CPWHLY00','CPASMA00','CCNUM00'),
           'mcs4_parent_cm_interview.dta'=c('DPWHLY00','DPASMA00','DCNUM00'),
           'mcs5_parent_cm_interview.dta'=c('EPWHLY00','EPASMA00','ECNUM00'),
           'mcs6_parent_cm_interview.dta'=c('FPASMA00','FCNUM00')),
  NextSteps=list('ns8_2015_main_interview.dta'=c('W8LOIL','W8LOLR0Z'),
                 'ns9_2022_main_interview.dta'=c('W9LOLRAD','W9DLOILCOND0A'))
)
sink(out)
for(cohort in names(target)) {
  cat('\n###',cohort,'\n')
  x <- readRDS(file.path(root,paste0(cohort,'_full.rds')))
  for(fn in names(target[[cohort]])) {
    ix <- which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn & x$source_inventory$format=='.dta')
    stopifnot(length(ix)==1L)
    d <- x$datasets[[ix]]
    cat('\nFILE',fn,'n=',nrow(d),'\n')
    for(v in target[[cohort]][[fn]]) {
      if(!v %in% names(d)) {cat('MISSING',v,'\n');next}
      col <- d[[v]]
      cat(v,' LABEL=',attr(col,'label'),'\n')
      print(attr(col,'labels'))
      print(table(as.character(col),useNA='always'))
    }
  }
  rm(x);gc()
}
sink()
cat(out,'\n')
