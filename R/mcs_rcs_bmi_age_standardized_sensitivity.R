input <- 'E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs'
d <- readRDS(file.path(input,'MCS_analysis_long.rds'))
d <- d[!is.na(d$wheeze12m) & is.finite(d$bmi),]
for(age in sort(unique(d$age))) {
 j <- d$age==age
 d$bmi_z[j] <- (d$bmi[j]-mean(d$bmi[j]))/sd(d$bmi[j])
}
d <- d[abs(d$bmi_z)<=2,]
d <- d[order(d$id,d$age),];d$id<-factor(d$id)
library(splines);library(geepack)
af<-function(z) ns(z,knots=7,Boundary.knots=c(3,11))
bf<-function(z) ns(z,knots=0,Boundary.knots=c(-2,2))
fit<-geeglm(wheeze12m~af(age)*bf(bmi_z),id=id,data=d,family=binomial,corstr='exchangeable')
co<-coef(fit);V<-vcov(fit);j<-grep(':',names(co),fixed=TRUE)
W<-as.numeric(t(co[j])%*%solve(V[j,j,drop=FALSE],co[j]))
z<-data.frame(n=nrow(d),children=length(unique(d$id)),Wald_chisq=W,df=length(j),p=pchisq(W,length(j),lower.tail=FALSE))
write.csv(z,file.path(input,'MCS_age_standardized_BMI_sensitivity.csv'),row.names=FALSE)
