root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/bmi_variable_audit.csv'
rows <- list(); k <- 0L
for (cohort in c('NCDS','BCS70','MCS','NextSteps')) {
  cat('Scanning',cohort,'\n'); flush.console()
  x <- readRDS(file.path(root,paste0(cohort,'_full.rds')))
  for (i in seq_along(x$datasets)) {
    d <- x$datasets[[i]]
    if (!is.data.frame(d)) next
    labs <- vapply(d,function(z) {a <- attr(z,'label'); if (is.null(a)) '' else as.character(a)[1]},'')
    hit <- which(grepl('(^|[^a-z])(bmi|body mass index|body-mass index|height|weight|obes|overweight)([^a-z]|$)',paste(names(d),labs),ignore.case=TRUE,perl=TRUE))
    if (!length(hit)) next
    fn <- basename(gsub('\\\\','/',x$source_inventory$source_entry[i]))
    for (j in hit) {
      z <- suppressWarnings(as.numeric(d[[j]]))
      k <- k+1L
      rows[[k]] <- data.frame(cohort=cohort,file=fn,variable=names(d)[j],label=labs[j],n=nrow(d),
         nonmissing=sum(!is.na(z)),positive=sum(is.finite(z)&z>0),
         min=if(any(is.finite(z)&z>0)) min(z[is.finite(z)&z>0]) else NA_real_,
         max=if(any(is.finite(z)&z>0)) max(z[is.finite(z)&z>0]) else NA_real_)
    }
  }
  rm(x);gc()
}
write.csv(do.call(rbind,rows),out,row.names=FALSE,fileEncoding='UTF-8')
cat(out,'rows',k,'\n')
