input <- 'E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs'
library(ggplot2)
d <- read.csv(file.path(input,'MCS_group_trajectory_observed.csv'))
s <- read.csv(file.path(input,'MCS_group_trajectory_class_summary.csv'))
labels <- c('C1 Early reported','C2 Later reported','C3 Rare reported','C4 Recurrent reported')
d$group <- factor(d$class,levels=1:4,labels=labels)
g <- ggplot(d,aes(age,wheeze12m,color=group))+geom_line(linewidth=1.1)+geom_point(size=2.5)+
 scale_x_continuous(breaks=c(3,5,7,11))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),labels=function(x) paste0(round(x*100),'%'))+
 scale_color_manual(values=c('#D27A41','#665B9B','#177A78','#BA4F58'))+
 labs(x='Nominal age (years)',y='Observed 12-month wheeze within assigned group',color=NULL,
 title='MCS four-class exploratory trajectory fit',subtitle='10,768 children with all four waves; class sizes: 1,689 / 482 / 7,368 / 1,229')+
 theme_classic(base_size=12)+theme(plot.title=element_text(face='bold',size=13),legend.position='bottom')
ggsave(file.path(input,'MCS_group_trajectories.png'),g,width=8,height=5,dpi=300,bg='white')
ggsave(file.path(input,'MCS_group_trajectories.pdf'),g,width=8,height=5,device=cairo_pdf,bg='white')
