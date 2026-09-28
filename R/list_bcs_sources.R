x<-readRDS('E:/公共数据库/Four_Cohorts_Full_RDS_20260925/BCS70_full.rds')
for(s in x$source_inventory$source_entry) cat(s,'\n')
