input <- 'E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs'
d <- readRDS(file.path(input,'MCS_analysis_long.rds'))
d <- d[!is.na(d$wheeze12m),c('id','age','wheeze12m')]
tab <- table(d$id)
d <- d[d$id %in% names(tab)[tab==4L],]
d <- d[order(d$id,d$age),]
d$id <- factor(d$id)
cat('complete subjects',nlevels(d$id),'rows',nrow(d),'\n')
library(flexmix)
set.seed(20260927)
fits <- list()
for(k in 2:4) {
 cat('fitting classes',k,'\n');flush.console()
 fits[[as.character(k)]] <- flexmix(cbind(wheeze12m,1-wheeze12m) ~ age + I(age^2) | id,data=d,k=k,model=FLXMRglm(family='binomial'),control=list(iter.max=150,tolerance=1e-6,verbose=0))
 saveRDS(fits[[as.character(k)]],file.path(input,paste0('MCS_group_trajectory_k',k,'.rds')),compress='xz')
 cat('BIC',BIC(fits[[as.character(k)]]),'\n');flush.console()
}
score <- data.frame(classes=2:4,BIC=vapply(fits,BIC,0.0),logLik=vapply(fits,function(z) as.numeric(logLik(z)),0.0))
write.csv(score,file.path(input,'MCS_group_trajectory_model_selection.csv'),row.names=FALSE)
best <- fits[[as.character(score$classes[which.min(score$BIC)])]]
post <- posterior(best)
subject <- unique(as.character(d$id))
if(nrow(post)==nrow(d)) post <- post[!duplicated(d$id),,drop=FALSE]
stopifnot(nrow(post)==length(subject))
cl <- max.col(post)
size <- as.data.frame(table(class=cl));names(size)[2]<-'children'
size$percent <- 100*size$children/sum(size$children)
size$mean_posterior <- vapply(seq_len(nrow(size)),function(j) mean(post[cl==j,j]),0.0)
write.csv(size,file.path(input,'MCS_group_trajectory_class_summary.csv'),row.names=FALSE)
d$class <- rep(cl,each=4L)
obs <- aggregate(wheeze12m~class+age,d,mean)
ns <- aggregate(wheeze12m~class+age,d,length);obs$n<-ns$wheeze12m
write.csv(obs,file.path(input,'MCS_group_trajectory_observed.csv'),row.names=FALSE)
saveRDS(list(fits=fits,best=best,selection=score),file.path(input,'MCS_group_trajectory_fits.rds'),compress='xz')
