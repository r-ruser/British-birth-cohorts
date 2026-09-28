x<-readRDS('E:/公共数据库/Four_Cohorts_Full_RDS_20260925/BCS70_full.rds')
fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
for(s in list(c('bcs7072a.dta','a0043b','a0255','a0278','a0014'),c('bcs1derived.dta','BD1PSOC'))) {
 d<-x$datasets[[which(fns==s[1])]];cat('\n',s[1],'\n')
 for(v in s[-1]) {z<-d[[v]];cat(v,'::',attr(z,'label'),'\n');print(attr(z,'labels'));if(length(unique(z))<40) print(table(z,useNA='always')) else print(summary(as.numeric(z)))}
}
