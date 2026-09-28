out<-'E:/cohort_rds_stage_20260925/formal_analysis_20260928'
library(ggplot2)
labels<-c(male='Male vs female',low_birthweight='Birthweight <2.5 kg',low_income='Lower family income',room_smoke='Smoking in child room',respondent_asthma='Parent respondent asthma',atopy3='Eczema/hay fever by 3',bmi3_z='BMI at 3, per SD',ever_wheeze3='Wheeze by 3',ever_asthma3='Asthma report by 3')

drawforest<-function(d,x,lo,hi,file,title,subtitle,xlim) {
 d$label<-factor(unname(labels[d$predictor]),levels=rev(unname(labels[d$predictor])))
 p<-ggplot(d,aes(x=.data[[x]],y=label))+geom_vline(xintercept=1,color='#808080',linetype='dashed')+
   geom_errorbarh(aes(xmin=.data[[lo]],xmax=.data[[hi]]),height=.18,color='#136F70',linewidth=.65)+
   geom_point(size=2.6,color='#136F70')+
   scale_x_log10(limits=xlim,breaks=c(.7,1,1.5,2,3),labels=c('0.7','1','1.5','2','3'))+
   labs(x='Adjusted ratio (95% CI)',y=NULL,title=title,subtitle=subtitle)+
   theme_classic(base_size=11)+theme(plot.title=element_text(face='bold',size=13),axis.text.y=element_text(size=10),plot.subtitle=element_text(size=10))
 ggsave(file.path(out,paste0(file,'.png')),p,width=7.2,height=5,dpi=300,bg='white')
 ggsave(file.path(out,paste0(file,'.pdf')),p,width=7.2,height=5,device=cairo_pdf,bg='white')
}
h<-read.csv(file.path(out,'MCS_weighted_interval_HR.csv'))
drawforest(h[h$model=='MCS_first_asthma_report',],'HR','CI_low','CI_high','MCS_first_asthma_HR_forest',
 'First subsequent asthma report','MCS age 3 to 14; weighted grouped-survival model',c(.65,3.5))
r<-read.csv(file.path(out,'MCS_age3_wheeze_to_age11_weighted_RR.csv'))
drawforest(r,'RR','CI_low','CI_high','MCS_adverse_prognosis_RR_forest',
 'Wheeze again at age 11','Among children with 12-month wheeze at age 3',c(.65,3.5))

a<-read.csv(file.path(out,'NCDS_BCS70_adult_interval_HR.csv'))
a<-a[a$method=='response_ipw' & a$predictor=='childhood_symptom',]
a$cohort<-factor(a$cohort,levels=c('NCDS','BCS70'))
p<-ggplot(a,aes(HR,cohort))+geom_vline(xintercept=1,color='#808080',linetype='dashed')+
 geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.12,linewidth=.7,color='#AE6840')+
 geom_point(size=3,color='#AE6840')+scale_x_log10(limits=c(.9,3.6),breaks=c(1,1.5,2,3))+
 labs(x='Adjusted interval HR (95% CI)',y=NULL,title='Childhood respiratory report and later adult report',
 subtitle='NCDS age 50 to 62; BCS70 age 38 to 46; response-weighted within cohort')+
 theme_classic(base_size=11)+theme(plot.title=element_text(face='bold',size=13))
ggsave(file.path(out,'NCDS_BCS70_childhood_adult_HR.png'),p,width=7.2,height=3.4,dpi=300,bg='white')
ggsave(file.path(out,'NCDS_BCS70_childhood_adult_HR.pdf'),p,width=7.2,height=3.4,device=cairo_pdf,bg='white')

c<-read.csv(file.path(out,'MCS_prognosis_calibration.csv'))
p<-ggplot(c,aes(predicted_oof,wheeze11))+geom_abline(slope=1,intercept=0,color='#8D8D8D',linetype='dashed')+
 geom_point(aes(size=n),color='#136F70',alpha=.85)+geom_line(color='#136F70',linewidth=.6)+
 scale_x_continuous(limits=c(0,.65),labels=function(x) paste0(round(100*x),'%'))+
 scale_y_continuous(limits=c(0,.65),labels=function(x) paste0(round(100*x),'%'))+
 labs(x='Predicted risk (out-of-fold)',y='Observed risk',title='Internal calibration of age-11 wheeze prediction',
 subtitle='Five-fold family-level cross-validation; points are predicted-risk deciles',size='Children')+
 theme_classic(base_size=11)+theme(plot.title=element_text(face='bold',size=13))
ggsave(file.path(out,'MCS_prognosis_calibration.png'),p,width=6.2,height=5.2,dpi=300,bg='white')
ggsave(file.path(out,'MCS_prognosis_calibration.pdf'),p,width=6.2,height=5.2,device=cairo_pdf,bg='white')
