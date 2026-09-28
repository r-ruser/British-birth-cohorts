out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
h<-read.csv(file.path(out,'MCS_weighted_interval_HR.csv'))
for(m in unique(h$model)) h$q_BH[h$model==m]<-p.adjust(h$p[h$model==m],method='BH')
write.csv(h,file.path(out,'MCS_weighted_interval_HR_with_FDR.csv'),row.names=FALSE)
r<-read.csv(file.path(out,'MCS_age3_wheeze_to_age11_weighted_RR.csv'))
r$q_BH<-p.adjust(r$p,method='BH')
write.csv(r,file.path(out,'MCS_age3_wheeze_to_age11_RR_with_FDR.csv'),row.names=FALSE)
