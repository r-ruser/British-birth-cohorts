x<-readRDS('E:/公共数据库/Four_Cohorts_Full_RDS_20260925/NCDS_full.rds')
ii<-which(basename(gsub('\\\\','/',x$source_inventory$source_entry))=='ncds_response.dta');stopifnot(length(ii)==1)
saveRDS(x$datasets[[ii]],'E:/cohort_rds_stage_20260925/life_course_v2_20260928/NCDS_response_source.rds')
