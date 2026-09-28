x<-readRDS('E:/公共数据库/Four_Cohorts_Full_RDS_20260925/BCS70_full.rds')
fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
for(fn in c('bcs7072a.dta','bcs7072b.dta','bcs1derived.dta')) {
 i<-which(fns==fn);d<-x$datasets[[i]];cat('\n',fn,'n',nrow(d),'\n')
 labs<-vapply(d,function(z) if(is.null(attr(z,'label'))) '' else as.character(attr(z,'label'))[1],'')
 j<-grep('bcsid|sex|birth.?weight|smok|social class|occupation|gestat|^.*weight.*baby',paste(names(d),labs),ignore.case=TRUE)
 for(k in j) cat(names(d)[k],'::',labs[k],'\n')
}
