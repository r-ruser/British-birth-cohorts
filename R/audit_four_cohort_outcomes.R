root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/four_cohort_outcome_variable_audit.csv'
terms <- c(
  asthma = 'asthma|asthmatic|wheez|bronchitis|respiratory',
  tuberculosis = 'tubercul|pulmonary tb|\btb\b|phthisis',
  lung_cancer = 'lung cancer|cancer of (the )?lung|bronchogenic|malignan.{0,30}lung|lung.{0,30}malignan'
)
hits <- list()
k <- 0L
for (cohort in c('BCS70','NCDS','MCS','NextSteps')) {
  cat('Scanning', cohort, '\n')
  bundle <- readRDS(file.path(root, paste0(cohort, '_full.rds')))
  for (j in seq_along(bundle$datasets)) {
    table <- bundle$datasets[[j]]
    if (!is.data.frame(table)) next
    var <- names(table)
    labs <- vapply(table, function(x) {
      z <- attr(x, 'label')
      if (is.null(z)) '' else as.character(z)[1]
    }, character(1))
    hay <- paste(var, labs)
    for (topic in names(terms)) {
      ix <- grep(terms[[topic]], hay, ignore.case=TRUE, perl=TRUE)
      if (!length(ix)) next
      k <- k + 1L
      hits[[k]] <- data.frame(
        cohort=cohort,
        topic=topic,
        source_entry=bundle$source_inventory$source_entry[j],
        source_key=bundle$source_inventory$key[j],
        rows=nrow(table),
        variable=var[ix],
        label=labs[ix],
        stringsAsFactors=FALSE
      )
    }
  }
  rm(bundle)
  gc()
}
result <- do.call(rbind, hits)
write.csv(result, out, row.names=FALSE, fileEncoding='UTF-8')
cat('Wrote', nrow(result), 'matches to', out, '\n')
print(aggregate(variable ~ cohort + topic, result, length))
