input <- 'E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs'
library(ggplot2)
p <- read.csv(file.path(input,'MCS_age_BMI_RCS_predictions.csv'))
w <- read.csv(file.path(input,'cohort_wave_prevalence.csv'))
w <- w[w$cohort=='MCS' & w$measure=='wheeze_last_12m',]
p$bmi_group <- factor(p$bmi,levels=sort(unique(p$bmi)),labels=c('BMI 10th percentile','BMI median','BMI 90th percentile'))
pal <- c('#177A78','#DA8246','#6C5B94')
g1 <- ggplot(p,aes(age,prob,color=bmi_group,fill=bmi_group))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.13,color=NA)+geom_line(linewidth=1.05)+
 geom_point(data=w,aes(x=age,y=prevalence),inherit.aes=FALSE,color='#202020',size=2.2)+
 scale_color_manual(values=pal)+scale_fill_manual(values=pal)+scale_x_continuous(breaks=c(3,5,7,11))+
 scale_y_continuous(labels=function(x) paste0(round(100*x),'%'),limits=c(0,.4))+
 labs(x='Nominal age (years)',y='Predicted 12-month wheeze',color=NULL,fill=NULL,
 title='Childhood wheeze across age and BMI',subtitle='MCS ages 3, 5, 7 and 11; dots show observed unweighted prevalence')+
 theme_classic(base_size=12)+theme(plot.title=element_text(face='bold',size=13),legend.position='bottom',legend.box='horizontal')
ggsave(file.path(input,'MCS_age_BMI_RCS_age_curves.png'),g1,width=8,height=5,dpi=300,bg='white')
ggsave(file.path(input,'MCS_age_BMI_RCS_age_curves.pdf'),g1,width=8,height=5,device=cairo_pdf,bg='white')

fit <- readRDS(file.path(input,'MCS_age_BMI_RCS_models.rds'))$interaction
library(splines)
ak <- 7; bk <- 16.6; bounds <- c(14.2,22.73)
af <- function(z) ns(z,knots=ak,Boundary.knots=c(3,11))
bf <- function(z) ns(z,knots=bk,Boundary.knots=bounds)
grid <- expand.grid(age=c(3,5,7,11),bmi=seq(bounds[1],bounds[2],length.out=100))
X <- model.matrix(delete.response(terms(fit)),grid)
eta <- as.vector(X%*%coef(fit));se <- sqrt(rowSums((X%*%vcov(fit))*X))
grid$prob <- plogis(eta);grid$lo <- plogis(eta-1.96*se);grid$hi <- plogis(eta+1.96*se)
grid$age<-factor(grid$age,levels=c(3,5,7,11),labels=paste0('Age ',c(3,5,7,11)))
g2 <- ggplot(grid,aes(bmi,prob,color=age,fill=age))+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.10,color=NA)+geom_line(linewidth=1)+
 scale_color_manual(values=c('#177A78','#DA8246','#6C5B94','#B84854'))+
 scale_fill_manual(values=c('#177A78','#DA8246','#6C5B94','#B84854'))+
 scale_y_continuous(labels=function(x) paste0(round(100*x),'%'))+
labs(x='BMI (kg/m^2)',y='Predicted 12-month wheeze',color=NULL,fill=NULL,
 title='BMI association differs by childhood age',subtitle='MCS; exploratory repeated-measures restricted cubic spline model')+
 theme_classic(base_size=12)+theme(plot.title=element_text(face='bold',size=13),legend.position='bottom')
ggsave(file.path(input,'MCS_age_BMI_RCS_BMI_curves.png'),g2,width=8,height=5,dpi=300,bg='white')
ggsave(file.path(input,'MCS_age_BMI_RCS_BMI_curves.pdf'),g2,width=8,height=5,device=cairo_pdf,bg='white')
