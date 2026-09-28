root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
x<-readRDS(file.path(root,'MCS_full.rds'))
fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry))
spec<-list('mcs1_hhgrid.dta'=c('AHCSEX00','ACNUM00'),
 'mcs1_cm_derived.dta'=c('ADBWGT00','ADGEST00','ACNUM00'),
 'mcs1_family_derived.dta'=c('AOECDUK0','ADWGT200','ADPTY200'),
 'mcs1_parent_interview.dta'=c('APNUM00','APSMKR00','APCIPR00','APSMCH00','APASTH00'),
 'mcs2_parent_cm_interview.dta'=c('BPECHA00','BPWHEE00','BPWHLY00','BPASMA00'),
 'mcs_longitudinal_cm_file.dta'=c('pttype2','weight2','bovwtgb','dovwtgb','enrespwt'))
for(fn in names(spec)) {
 i<-which(fns==fn);d<-x$datasets[[i]];cat('\n',fn,'n=',nrow(d),'\n')
 for(v in spec[[fn]]) {
  if(!v %in% names(d)) {cat('MISSING',v,'\n');next}
  z<-d[[v]];cat(v,'::',attr(z,'label'),'\n')
  print(attr(z,'labels'))
  if(length(unique(z))<30) print(table(z,useNA='always')) else print(summary(as.numeric(z)))
 }
}
