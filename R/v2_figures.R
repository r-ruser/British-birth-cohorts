library(ggplot2);library(patchwork);library(survey)
out<-'E:/cohort_rds_stage_20260925/life_course_v2_20260928';dir.create(file.path(out,'figures'),showWarnings=FALSE)
theme_set(theme_classic(base_size=10,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',plot.title=element_text(face='bold',size=12)))
cols<-c(MCS='#557D9E',NCDS='#B2816B',BCS70='#73927D')
savep<-function(p,name,w=9,h=5){ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=320,bg='white');ggsave(file.path(out,'figures',paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf);ggsave(file.path(out,'figures',paste0(name,'.svg')),p,width=w,height=h,device=svglite::svglite)}
r<-read.csv(file.path(out,'model_RR_results.csv'));i<-read.csv(file.path(out,'interaction_tests.csv'));f<-read.csv(file.path(out,'attrition_flow_weights.csv'));s<-read.csv(file.path(out,'sex_age_adjusted_RR.csv'))
labs<-c(male='Male vs female',low_birthweight='Birthweight <2,500 g',ses='Early socioeconomic disadvantage',pregnancy_continued_smoke='Prenatal smoking proxy',respondent_asthma='Parent asthma',atopy3='Atopy at age 3',bmi3_z='BMI at age 3 (per SD)',ever_asthma3='Ever-asthma report by age 3',childhood_symptom='Childhood respiratory symptoms')
forest<-function(z){z$predictor<-factor(labs[z$term],levels=rev(unique(labs[z$term])));ggplot(z,aes(estimate,predictor,color=cohort))+geom_vline(xintercept=1,lty=2,color='grey65')+geom_errorbar(aes(xmin=lower,xmax=upper),width=.15,orientation='y')+geom_point(size=2)+scale_x_log10()+scale_color_manual(values=cols)+labs(x='Adjusted RR (95% CI)',y=NULL,color=NULL)}
z<-r[r$method=='IPCW' & r$endpoint %in% c('transition_0_to_1','transition_1_to_1') & r$term %in% names(labs),];z$transition<-ifelse(z$endpoint=='transition_0_to_1','No recent wheeze to recent wheeze','Recent wheeze to recent wheeze')
savep(forest(z)+facet_wrap(~transition,nrow=1)+labs(title='Early susceptibility differs for symptom appearance and persistence'),'01_MCS_transitions',11,5)
z<-r[r$method=='IPCW' & r$endpoint=='adult_repeated_burden' & r$term %in% names(labs),]
savep(forest(z)+facet_wrap(~cohort,nrow=1)+labs(title='Childhood symptoms remain associated with adult reported respiratory burden'),'02_adult_GEE',10,4.4)
z<-r[r$method=='response_IPW' & r$endpoint=='shared_age_10_11' & r$term %in% names(labs),]
savep(forest(z)+facet_wrap(~cohort,nrow=1)+labs(title='Age 10/11 replication across birth cohorts'),'03_common_age_replication',10,4.2)
ss<-s[s$method=='IPCW',];p<-ggplot(ss,aes(age,RR,color=cohort))+geom_hline(yintercept=1,lty=2,color='grey60')+geom_ribbon(aes(ymin=lower,ymax=upper,fill=cohort),alpha=.16,color=NA)+geom_line()+geom_point(size=2)+facet_wrap(~cohort,scales='free_x')+scale_color_manual(values=cols)+scale_fill_manual(values=cols)+labs(x='Survey age (years)',y='Adjusted male / female RR (95% CI)',color=NULL,fill=NULL,title='Sex differences vary across observed ages')
savep(p,'04_sex_age',10,4.4)
ff<-reshape(f[c('cohort','age','baseline_complete','monotone_observed','all_available')],varying=list(c('monotone_observed','all_available')),v.names='n',timevar='sample',times=c('Monotone follow-up','All observed visits'),direction='long');ff$percent<-100*ff$n/ff$baseline_complete
p<-ggplot(ff,aes(age,percent,linetype=sample,color=cohort))+geom_line()+geom_point()+facet_wrap(~cohort,scales='free_x')+scale_color_manual(values=cols)+labs(x='Survey age (years)',y='Observed / baseline complete cases (%)',linetype=NULL,color=NULL,title='Attrition IPCW targets a defined monotone follow-up population')
savep(p,'05_attrition',10,4.4)
ll<-readRDS(file.path(out,'analysis_longitudinal_frames.rds'));z<-ll$MCS$monotone;z<-z[z$age>3,];z$interval<-factor(z$age,levels=c(5,7,11),labels=c('3-5 years','5-7 years','7-11 years'));z$origin<-factor(z$previous,labels=c('Previously no recent wheeze','Previously recent wheeze'));des<-svydesign(ids=~family,strata=~design_stratum,weights=~w,data=z,nest=TRUE)
tt<-svyby(~y,~interval+origin,des,svymean,vartype='ci');write.csv(tt,file.path(out,'MCS_transition_absolute_risks.csv'),row.names=FALSE)
p<-ggplot(tt,aes(interval,y,color=origin,group=origin))+geom_line()+geom_point(size=2)+geom_errorbar(aes(ymin=ci_l,ymax=ci_u),width=.1)+scale_y_continuous(labels=scales::label_percent())+scale_color_manual(values=c('#557D9E','#B2816B'))+labs(x=NULL,y='Next-wave recent wheeze (%)',color=NULL,title='Prior recent wheeze identifies a high reporting-risk group')
savep(p,'06_MCS_absolute_transitions',8,4.3)
ww<-do.call(rbind,lapply(names(ll),function(c){z<-ll[[c]]$monotone;z<-z[z$age==max(z$age),];data.frame(cohort=c,IPCW=z$ipcw_raw)}))
p<-ggplot(ww,aes(IPCW,fill=cohort))+geom_histogram(bins=45,color='white',linewidth=.15)+facet_wrap(~cohort,scales='free')+scale_fill_manual(values=cols)+labs(x='Untrimmed stabilized IPCW at final wave',y='Participants',fill=NULL,title='Estimated dropout weights have limited extreme values')
savep(p,'07_weight_diagnostics',10,4)
writeLines(c('Figure contract: quantitative grids; evidence chain susceptibility to conditional transitions to adult burden to birth-cohort replication.', 'R-only rendering; white background; Arial; 320 dpi PNG plus editable SVG and PDF.', 'Forest CIs: survey or cluster-robust Wald 95%; sex-age ribbons: adjusted contrast CIs; transition absolute-risk CIs: survey linearization.', 'Attrition/weight figures are descriptive diagnostics; primary weights are trimmed at wave-specific 1st/99th percentiles.'),file.path(out,'figure_contract.txt'))
cat('Seven figures complete\n')
