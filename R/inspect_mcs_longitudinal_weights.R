root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
x<-readRDS(file.path(root,'MCS_full.rds'))
fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
for(fn in c('mcs_longitudinal_cm_file.dta','mcs_core.dta')) {
 i<-which(fns==fn);d<-x$datasets[[i]]
 cat('\n',fn,'n=',nrow(d),'\n')
 cat('WEIGHT COLUMNS',paste(grep('ovwt|nresp|weight',names(d),value=TRUE,ignore.case=TRUE),collapse='|'),'\n')
 labs<-vapply(d,function(z) if(is.null(attr(z,'label'))) '' else as.character(attr(z,'label'))[1],'')
 j<-grep('mcsid|cnum|weight|wgt|wtgb|nresp|strat|sex|birth weight|income quint',paste(names(d),labs),ignore.case=TRUE)
 for(k in j) if(grepl('mcsid|cnum|weight.*use|wgt|wtgb|nresp|strat|^.*sex|birth weight|income quint',paste(names(d)[k],labs[k]),ignore.case=TRUE)) {
  z<-d[[k]];cat(names(d)[k],'::',labs[k],':: valid',sum(!is.na(z)&suppressWarnings(as.numeric(z))>=0),'\n')
 }
}
