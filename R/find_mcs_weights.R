root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
x <- readRDS(file.path(root,'MCS_full.rds'))
for(i in seq_along(x$datasets)) {
 d<-x$datasets[[i]];if(!is.data.frame(d)) next
 labs<-vapply(d,function(z) if(is.null(attr(z,'label'))) '' else as.character(attr(z,'label'))[1],'')
 j<-which(grepl('weight for use|survey weight|sampling weight|non.response weight|longitudinal weight|^[A-Za-z]{1,4}WGT[0-9A-Z_]*$',paste(names(d),labs),ignore.case=TRUE))
 if(length(j)) {cat(basename(gsub('\\\\','/',x$source_inventory$source_entry[i])),'\n'); for(k in j) cat(' ',names(d)[k],labs[k],'\n')}
}
