out <- 'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
x <- readRDS(file.path(out, 'MCS_asthma_first_report_intervals.rds'))
covars <- c('male','low_birthweight','low_income','room_smoke',
            'respondent_asthma','atopy3','bmi3_z','ever_wheeze3')
x <- x[complete.cases(x[c('event','weight','family','stratum',covars)]), ]
groups <- c('male','low_birthweight','low_income','respondent_asthma','atopy3','ever_wheeze3')
ages <- c('3-5','5-7','7-11','11-14')
z <- list(); k <- 0L
for (v in groups) for (g in 0:1) {
  a <- x[x[[v]] == g, ]
  s <- 1
  for (interval in ages) {
    b <- a[a$interval == interval, ]
    if (nrow(b) == 0) next
    p <- weighted.mean(b$event, b$weight)
    s <- s * (1-p)
    k <- k+1L
    z[[k]] <- data.frame(group=v,level=g,interval=interval,
                         at_risk=nrow(b),events=sum(b$event),
                         weighted_interval_risk=p,
                         cumulative_report_risk=1-s)
  }
}
write.csv(do.call(rbind,z),file.path(out,'MCS_group_absolute_report_risk.csv'),row.names=FALSE)
