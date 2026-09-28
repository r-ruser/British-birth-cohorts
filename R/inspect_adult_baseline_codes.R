root<-'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
for(cohort in c('NCDS','BCS70')) {
 x<-readRDS(file.path(root,paste0(cohort,'_full.rds')))
 fn<-if(cohort=='NCDS') 'ncds0123.dta' else 'sn3723.dta'
 i<-which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn)
 d<-x$datasets[[i]]
 vs<-if(cohort=='NCDS') c('n622','n503','n639','n492','n646','n259','n260') else c('sex10','c3_4','b7_1','b7_2')
 cat('\n',cohort,'\n')
 for(v in vs) {z<-d[[v]];cat(v,'::',attr(z,'label'),'\n');print(attr(z,'labels'));print(sort(table(z,useNA='always'),decreasing=TRUE))}
}
