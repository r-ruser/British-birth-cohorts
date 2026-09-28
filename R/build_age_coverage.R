input <- 'E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs'
w <- read.csv(file.path(input,'cohort_wave_prevalence.csv'))
groups <- data.frame(group=c('newborn','infant','young child','older child','early adolescent','older adolescent','adolescent','youth','young people','adult 25-34','adult 35-44','adult 45-54','adult 55-64'),
 min_age=c(0,29/365,1,5,10,15,10,15,10,25,35,45,55),
 max_age=c(28/365,1,5,10,15,20,20,25,25,35,45,55,65),
 overlaps=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE))
z <- list()
for(cohort in unique(w$cohort)) for(i in seq_len(nrow(groups))) {
 a<-unique(w$age[w$cohort==cohort & w$age>=groups$min_age[i] & w$age<groups$max_age[i]])
 z[[length(z)+1L]]<-data.frame(cohort,group=groups$group[i],age_range=paste0(groups$min_age[i],'-',groups$max_age[i]),overlapping_descriptor=groups$overlaps[i],observed_ages=if(length(a)) paste(sort(a),collapse='|') else '',n_distinct_nominal_ages=length(a))
}
write.csv(do.call(rbind,z),file.path(input,'age_group_coverage.csv'),row.names=FALSE,fileEncoding='UTF-8')
