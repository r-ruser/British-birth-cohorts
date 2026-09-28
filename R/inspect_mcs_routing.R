root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
x <- readRDS(file.path(root,'MCS_full.rds'))
for (p in c('B','C','D','E')) {
 fn <- switch(p,B='mcs2_parent_cm_interview.dta',C='mcs3_parent_cm_interview.dta',D='mcs4_parent_cm_interview.dta',E='mcs5_parent_cm_interview.dta')
 ix <- which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn)
 d<-x$datasets[[ix]]
 cat('\n',p,'\n');print(attr(d[[paste0(p,'PWHEE00')]],'labels'))
 print(table(d[[paste0(p,'PWHEE00')]],d[[paste0(p,'PWHLY00')]],useNA='always'))
 id <- paste(d$MCSID,d[[paste0(p,'CNUM00')]],sep='-')
 valid <- d[[paste0(p,'PWHEE00')]] %in% c(1,2)
 cat('valid ever-wheeze duplicates',sum(duplicated(id[valid])),' unique valid',length(unique(id[valid])),'\n')
}
