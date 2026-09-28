root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pattern <- 'smok|cigarette|breast|b.?feed|alcohol|drink|parent.*read|read.*child|bed.?time|disciplin|parent.*activity|mother.*behav|father.*behav'
rows <- list(); k <- 0L
for(cohort in c('MCS','NCDS','BCS70','NextSteps')) {
 cat('audit ',cohort,'\n');flush.console()
 x <- readRDS(file.path(root,paste0(cohort,'_full.rds')))
 for(i in seq_along(x$datasets)) {
  d <- x$datasets[[i]];if(!is.data.frame(d)) next
  fn <- basename(gsub('\\\\','/',x$source_inventory$source_entry[i]))
  include <- switch(cohort,MCS=grepl('mcs1_|mcs2_',fn),NCDS=grepl('ncds0123',fn),
                     BCS70=grepl('bcs7072a|bcs1derived|sn3723|bcs70_5',fn),NextSteps=grepl('wave1|wave2|lsype.*1|lsype.*2|lsype1|lsype2',fn,ignore.case=TRUE))
  if(!include) next
  labs<-vapply(d,function(z){a<-attr(z,'label');if(is.null(a)) '' else as.character(a)[1]},'')
  for(j in which(grepl(pattern,paste(names(d),labs),ignore.case=TRUE,perl=TRUE))) {
   z<-d[[j]];k<-k+1L
   vl<-attr(z,'labels');tb<-sort(table(as.character(z),useNA='ifany'),decreasing=TRUE)
   rows[[k]]<-data.frame(cohort,file=fn,variable=names(d)[j],label=labs[j],n=nrow(d),
    values=if(is.null(vl)) '' else paste(names(vl),vl,sep=':',collapse='|'),
    frequencies=paste(names(head(tb,20)),head(tb,20),sep=':',collapse='|'))
  }
 }
 rm(x);gc()
}
write.csv(do.call(rbind,rows),file.path(out,'parent_behavior_candidate_audit.csv'),row.names=FALSE,fileEncoding='UTF-8')
cat('rows ',k,'\n')
