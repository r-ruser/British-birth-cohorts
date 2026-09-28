out <- 'E:/cohort_rds_stage_20260925/parent_behavior_20260928'
library(ggplot2)
r<-read.csv(file.path(out,'parent_behavior_HR_RR.csv'))
nm<-c(pregnancy_continued_smoke='Pregnancy smoking\n(continuation proxy)',mother_smoke='Mother smoking',father_smoke='Father smoking',any_parent_smoke='Any parent smoking',breastfed='Any breastfeeding')
r$behavior<-factor(nm[r$exposure],levels=rev(nm))
r$adjustment<-factor(r$adjustment,levels=c('basic','extended'),labels=c('Basic adjustment','Extended adjustment'))
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(legend.position='bottom',legend.title=element_blank(),strip.background=element_blank(),strip.text=element_text(face='bold',size=9),plot.title=element_text(face='bold'),axis.text=element_text(colour='black'),panel.grid=element_blank()))
forest<-function(z,title,subtitle,xlab) {
 ggplot(z,aes(estimate,behavior,colour=adjustment))+geom_vline(xintercept=1,linetype=2,colour='grey55',linewidth=.35)+
 geom_errorbar(aes(xmin=CI_low,xmax=CI_high),orientation='y',width=.25,position=position_dodge(width=.5),linewidth=.45)+
 geom_point(position=position_dodge(width=.5),size=1.8)+
 scale_x_log10(breaks=c(.5,.75,1,1.5,2),labels=c('0.5','0.75','1','1.5','2'))+
 scale_colour_manual(values=c('#247B83','#8B8B8B'))+
 labs(x=xlab,y=NULL,title=title,subtitle=subtitle)
}
savefig<-function(p,name,w,h) {
 ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white')
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)
 ggsave(file.path(out,paste0(name,'.svg')),p,width=w,height=h,device=svglite::svglite)
}
z<-r[r$scale=='HR' & r$endpoint!='first_wheeze_report',]
z$panel<-factor(z$cohort,levels=c('MCS','NCDS','BCS70'),labels=c('MCS: age 3 to 14\nFirst asthma report','NCDS: age 50 to 62\nFirst adult re-report','BCS70: age 38 to 46\nFirst adult re-report'))
p<-forest(z,'Parental behaviors and later respiratory reports','Cohort-specific estimates; exposure timing and outcome definitions differ','Adjusted interval HR (95% CI)')+facet_grid(.~panel)
savefig(p,'parent_behavior_HR_by_cohort',10.2,4.8)
z<-r[r$scale=='RR',]
p<-forest(z,'Wheeze again at age 11','MCS children with current wheeze at age 3; smoking estimates are conditional associations','Adjusted RR (95% CI)')
savefig(p,'parent_behavior_MCS_prognosis_RR',6.8,4.7)
