out<-'E:/cohort_rds_stage_20260925/life_course_v3_20260928'
r<-read.csv(file.path(out,'MI_IPCW_core_results.csv'));lambda<-if('lambda'%in%names(r))r$lambda else r$FMI;r$lambda<-lambda;r$FMI<-lambda+(1-lambda)*2/(r$df+3)
write.csv(r,file.path(out,'MI_IPCW_core_results.csv'),row.names=FALSE)
t<-readLines('v3_mi.R');t<-gsub('FMI=\\(1\\+1/m\\)\\*between/total','lambda=(1+1/m)*between/total,FMI=((1+1/m)*between/total)+(1-(1+1/m)*between/total)*2/(df+3)',t);writeLines(t,'v3_mi.R')
cat('Rubin lambda and FMI separately reported\n')
