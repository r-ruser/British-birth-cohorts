root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/formal_covariate_audit.csv'
rows <- list(); k <- 0L
pattern <- 'sex|gender|birth.?weight|gestat|prematur|maternal.*smok|smok.*preg|smok.*mother|parent.*asth|mother.*asth|father.*asth|family.*asth|social class|socio.?economic|depriv|income|education|occupat|ethnic|survey weight|longitudinal weight|non.?response weight|design weight|strat|cluster|sampling weight'
for(cohort in c('MCS','NCDS','BCS70','NextSteps')) {
 cat('audit',cohort,'\n');flush.console()
 x<-readRDS(file.path(root,paste0(cohort,'_full.rds')))
 inv<-x$source_inventory
 for(i in seq_along(x$datasets)) {
  d<-x$datasets[[i]];if(!is.data.frame(d)) next
  fn<-basename(gsub('\\\\','/',inv$source_entry[i]))
  # Prioritize origin/baseline, outcomes and derived or weight files.
  include<-switch(cohort,
   MCS=grepl('mcs1_|mcs2_|mcs5_|mcs6_|mcs_core|weight|response',fn,ignore.case=TRUE),
   NCDS=grepl('ncds0123|ncds_2008|ncds_2013|ncds10|weight',fn,ignore.case=TRUE),
   BCS70=grepl('sn3723|bcs7016x|bcs_2008|bcs70_2012|bcs_age46|weight',fn,ignore.case=TRUE),
   NextSteps=grepl('ns8_2015|ns9_2022|weight',fn,ignore.case=TRUE))
  if(!include) next
  labs<-vapply(d,function(z) {a<-attr(z,'label');if(is.null(a)) '' else as.character(a)[1]},'')
  hit<-which(grepl(pattern,paste(names(d),labs),ignore.case=TRUE,perl=TRUE))
  if(!length(hit)) next
  for(j in hit) {
   k<-k+1L;z<-d[[j]];vals<-suppressWarnings(as.numeric(z)); pos<-vals[is.finite(vals)&vals>=0]
   rows[[k]]<-data.frame(cohort,file=fn,variable=names(d)[j],label=labs[j],n=nrow(d),n_nonnegative=length(pos),
      n_unique=length(unique(pos)),value_labels=if(is.null(attr(z,'labels'))) '' else paste(names(attr(z,'labels')),attr(z,'labels'),sep=':',collapse='|'))
  }
 }
 rm(x);gc()
}
write.csv(do.call(rbind,rows),out,row.names=FALSE,fileEncoding='UTF-8')
cat('wrote',out,'rows',k,'\n')
