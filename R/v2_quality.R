out<-'E:/cohort_rds_stage_20260925/life_course_v2_20260928'
x<-readRDS(file.path(out,'selected_source_tables.rds'))
sink(file.path(out,'routing_crosschecks.txt'))
for(spec in list(c('BCS70:bcs96x.dta','b960448','b960514'),c('NCDS:ncds4.dta','n5770','n5771'),c('NCDS:ncds0123.dta','n1305','n1306'),c('BCS70:sn3723.dta','b7_1','b7_20'))){d<-x[[spec[1]]];cat('\n',paste(spec,collapse=' / '),'\n');print(table(d[[spec[2]]],d[[spec[3]]],useNA='always'))}
sink()
r<-read.csv(file.path(out,'model_RR_results.csv'));f<-read.csv(file.path(out,'attrition_flow_weights.csv'))
stopifnot(all(is.finite(r$estimate)),all(r$lower<=r$estimate),all(r$estimate<=r$upper),all(f$events<=f$monotone_observed),!any(f$all_available<f$monotone_observed))
stopifnot(f$all_available[f$cohort=='BCS70'&f$age==26]>6000)
writeLines('PASS: finite estimates, valid CI ordering, outcome counts, monotone risksets, age26 routing denominator.',file.path(out,'quality_checks.txt'))
# Audit individual vital status in the actual bundle, without assuming a death registry.
root<-'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
sink(file.path(out,'vital_status_audit.txt'))
for(c in c('NCDS','BCS70')){
 cat('\n',c,'\n');b<-readRDS(file.path(root,paste0(c,'_full.rds')));fns<-basename(gsub('\\\\','/',b$source_inventory$source_entry))
 for(i in seq_along(b$datasets)){
  d<-b$datasets[[i]];if(!is.data.frame(d))next
  labs<-vapply(d,function(v){l<-attr(v,'label');if(is.null(l))''else as.character(l)[1]},'')
  ii<-grep('cohort member.*(died|dead|death)|death.*cohort member|date.*death|^CM.*(died|dead)|vital status|deceased',labs,ignore.case=TRUE)
  if(grepl('response',fns[i]))ii<-unique(c(ii,seq_along(labs)))
  for(j in ii)cat(fns[i],names(d)[j],labs[j],paste(names(attr(d[[j]],'labels')),attr(d[[j]],'labels'),collapse=' | '),'\n')
 }
 rm(b);gc()
}
sink()
cat('QA complete\n')
