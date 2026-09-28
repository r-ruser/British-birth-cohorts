options(stringsAsFactors=FALSE)
root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
out <- 'E:/cohort_rds_stage_20260925/four_cohort_trajectory_rcs'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
getd <- function(x,fn) {
 i <- which(basename(gsub('\\\\','/',x$source_inventory$source_entry))==fn & x$source_inventory$format=='.dta')
 stopifnot(length(i)==1L)
 x$datasets[[i]]
}
num <- function(x) suppressWarnings(as.numeric(x))
yn <- function(x,yes=1,no=2) ifelse(num(x)==yes,1L,ifelse(num(x)==no,0L,NA_integer_))
save_csv <- function(x,fn) write.csv(x,file.path(out,fn),row.names=FALSE,fileEncoding='UTF-8')
summ <- list(); trans <- list(); cohort_data <- list()
add_summary <- function(cohort,age,stage,measure,y,bmi=rep(NA_real_,length(y))) {
 ok<-!is.na(y)
 summ[[length(summ)+1L]] <<- data.frame(cohort,age,stage,measure,n_eligible=length(y),n_valid=sum(ok),n_yes=sum(y==1,na.rm=TRUE),prevalence=if(sum(ok)) mean(y[ok]) else NA_real_,n_bmi_valid=sum(ok & is.finite(bmi) & bmi>=8 & bmi<=60))
}
add_transition <- function(cohort,measure,d,ages) {
 for(i in seq_len(length(ages)-1L)) {
  a<-d[[paste0('y',ages[i])]]; b<-d[[paste0('y',ages[i+1L])]]
  z<-!is.na(a)&!is.na(b)
  for(v in 0:1) for(w in 0:1) trans[[length(trans)+1L]] <<- data.frame(cohort,measure,age_from=ages[i],age_to=ages[i+1L],from=v,to=w,n=sum(z & a==v & b==w))
 }
}
stage <- function(age) ifelse(age<1,'infant',ifelse(age<5,'young_child_1_4',ifelse(age<10,'older_child_5_9',ifelse(age<15,'early_adolescent_10_14',ifelse(age<20,'late_adolescent_15_19',ifelse(age<25,'adult_20_24',ifelse(age<35,'adult_25_34',ifelse(age<45,'adult_35_44',ifelse(age<55,'adult_45_54','adult_55_64')))))))))

cat('MCS\n'); x<-readRDS(file.path(root,'MCS_full.rds'))
mm<-list()
spec<-data.frame(age=c(3,5,7,11,14),p=c('B','C','D','E','F'),
 parent=c('mcs2_parent_cm_interview.dta','mcs3_parent_cm_interview.dta','mcs4_parent_cm_interview.dta','mcs5_parent_cm_interview.dta','mcs6_parent_cm_interview.dta'),
 bmi_file=c('mcs2_cm_interview.dta','mcs3_cm_interview.dta','mcs4_cm_interview.dta','mcs5_cm_interview.dta','mcs6_cm_derived.dta'),
 bmi_var=c('BCBMIN00','BMIN3','DCBMIN4','EBMIN5','FCBMIN6'))
for(i in seq_len(nrow(spec))) {
 s<-spec[i,]; p<-s$p; age<-s$age
 a<-getd(x,s$parent); b<-getd(x,s$bmi_file)
 idcol<-paste0(p,'CNUM00')
 aa<-data.frame(id=paste(a$MCSID,num(a[[idcol]]),sep=':'),ever_asthma=yn(a[[paste0(p,'PASMA00')]]))
 if(age<=11) {
  ev<-num(a[[paste0(p,'PWHEE00')]]); cur<-num(a[[paste0(p,'PWHLY00')]])
  aa$wheeze12m<-ifelse(ev==2,0L,ifelse(ev==1 & cur==1,1L,ifelse(ev==1 & cur==2,0L,NA_integer_)))
 }
 # Parent table may contain additional records for other respondents. Keep the sole valid response per child.
 aa$valid_score<-as.integer(!is.na(aa$ever_asthma))+if('wheeze12m'%in%names(aa)) as.integer(!is.na(aa$wheeze12m)) else 0L
 aa<-aa[order(aa$id,-aa$valid_score),];aa<-aa[!duplicated(aa$id),];aa$valid_score<-NULL
 bb<-data.frame(id=paste(b$MCSID,num(b[[idcol]]),sep=':'),bmi=num(b[[s$bmi_var]]))
 bb<-bb[!duplicated(bb$id),]
 d<-merge(aa,bb,by='id',all.x=TRUE,sort=FALSE)
 d$bmi[d$bmi<8 | d$bmi>60]<-NA_real_;d$age<-age
 add_summary('MCS',age,stage(age),'reported_ever_asthma',d$ever_asthma,d$bmi)
 if(age<=11) add_summary('MCS',age,stage(age),'wheeze_last_12m',d$wheeze12m,d$bmi)
 mm[[as.character(age)]]<-d
 cat(' age',age,'rows',nrow(d),'bmi',sum(!is.na(d$bmi)),'\n')
}
for(measure in c('ever_asthma','wheeze12m')) {
 ages<-if(measure=='ever_asthma') c(3,5,7,11,14) else c(3,5,7,11)
 wide<-data.frame(id=unique(unlist(lapply(mm[as.character(ages)],function(d)d$id))))
 for(age in ages) {d<-mm[[as.character(age)]]; tmp<-d[c('id',measure)];names(tmp)[2]<-paste0('y',age);wide<-merge(wide,tmp,by='id',all.x=TRUE,sort=FALSE)}
 add_transition('MCS',if(measure=='ever_asthma') 'reported_ever_asthma' else measure,wide,ages)
 complete<-complete.cases(wide[paste0('y',ages)])
 pat<-apply(wide[complete,paste0('y',ages),drop=FALSE],1,paste0,collapse='')
 tab<-sort(table(pat),decreasing=TRUE)
 save_csv(data.frame(pattern=names(tab),n=as.integer(tab),complete_case_n=sum(complete)),paste0('MCS_',measure,'_patterns.csv'))
 if(measure=='wheeze12m') {
  z<-wide[complete,]; y<-z[paste0('y',ages)]
  g<-ifelse(rowSums(y)==0,'never',ifelse(rowSums(y)==length(ages),'persistent_all',ifelse(y[[1]]==1 & y[[length(ages)]]==0,'early_then_absent',ifelse(y[[1]]==0 & y[[length(ages)]]==1,'later_present','intermittent'))))
  save_csv(as.data.frame(table(group=g)), 'MCS_wheeze_descriptive_groups.csv')
 }
}
mcs_long<-do.call(rbind,lapply(mm[c('3','5','7','11')],function(d)d[c('id','age','wheeze12m','ever_asthma','bmi')]))
saveRDS(mcs_long,file.path(out,'MCS_analysis_long.rds'),compress='xz')
rm(x,mm);gc()

cat('NCDS\n');x<-readRDS(file.path(root,'NCDS_full.rds'))
a<-getd(x,'ncds0123.dta'); ids<-as.character(a$ncdsid)
asth7<-yn(a$n259,2,3); bron7<-yn(a$n260,2,3)
any7<-ifelse(asth7==1 | bron7==1,1L,ifelse(asth7==0 & bron7==0,0L,NA_integer_))
v11<-num(a$n1305); v16<-num(a$n2617)
nch<-data.frame(id=ids,y7=any7,y11=ifelse(v11 %in% c(1,2),1L,ifelse(v11==3,0L,NA_integer_)),y16=yn(v16,1,2))
for(age in c(7,11,16)) add_summary('NCDS',age,stage(age),'reported_ever_asthma_or_wheezy_bronchitis',nch[[paste0('y',age)]])
add_transition('NCDS','reported_ever_asthma_or_wheezy_bronchitis',nch,c(7,11,16))
adult_spec<-data.frame(age=c(50,55,62),file=c('ncds_2008_followup.dta','ncds_2013_flatfile.dta','ncds10_age62_main_interview.dta'),v=c('N8KHPB01','N9KHPB01','n10khlprb01'),no=c(0,2,0),bmi_file=c('','ncds_2013_derived.dta','ncds10_age62_main_interview.dta'),bmi_var=c('','ND9BMI','nd10bmi_rec'))
nad<-data.frame(id=unique(ids))
for(i in seq_len(nrow(adult_spec))) {
 s<-adult_spec[i,];d<-getd(x,s$file);y<-yn(d[[s$v]],1,s$no); id<-as.character(d[[grep('^ncdsid$',names(d),ignore.case=TRUE,value=TRUE)[1]]]);
 bmi<-rep(NA_real_,length(y))
 if(nchar(s$bmi_var)) {e<-getd(x,s$bmi_file);j<-match(id,as.character(e[[grep('^ncdsid$',names(e),ignore.case=TRUE,value=TRUE)[1]]]));bmi<-num(e[[s$bmi_var]][j])}
 add_summary('NCDS',s$age,stage(s$age),'adult_asthma_or_wheezy_bronchitis',y,bmi)
 q<-data.frame(id=id,y=y);names(q)[2]<-paste0('y',s$age);nad<-merge(nad,q,by='id',all.x=TRUE,sort=FALSE)
}
add_transition('NCDS','adult_asthma_or_wheezy_bronchitis',nad,c(50,55,62))
rm(x,a,nch,nad);gc()

cat('BCS70\n');x<-readRDS(file.path(root,'BCS70_full.rds'))
a<-getd(x,'sn3723.dta');b<-getd(x,'bcs7016x.dta')
bc<-data.frame(id=as.character(a$bcsid),y10=yn(a$b7_1,1,2))
q<-data.frame(id=as.character(b$bcsid),y16=yn(b$ob9_1,1,2));bc<-merge(bc,q,by='id',all.x=TRUE,sort=FALSE)
add_summary('BCS70',10,stage(10),'reported_wheeze_attacks',bc$y10)
add_summary('BCS70',16,stage(16),'reported_wheeze_attacks',yn(b$ob9_1,1,2))
add_transition('BCS70','reported_wheeze_attacks',bc,c(10,16))
add_summary('BCS70',16,stage(16),'asthma_last_12m',yn(b$ha4_6,1,2))
bad<-data.frame(id=unique(as.character(a$bcsid)))
bs<-data.frame(age=c(38,42,46),file=c('bcs_2008_followup.dta','bcs70_2012_flatfile.dta','bcs_age46_main.dta'),v=c('b8khpb01','B9KHPB02','B10KHPB01'),bmi_file=c('','bcs70_2012_derived.dta','bcs_age46_main.dta'),bmi_var=c('','BD9BMI','BD10BMI'))
for(i in seq_len(nrow(bs))) {
 s<-bs[i,];d<-getd(x,s$file);id<-as.character(d[[grep('^bcsid$',names(d),ignore.case=TRUE,value=TRUE)[1]]]);y<-yn(d[[s$v]],1,0)
 bmi<-rep(NA_real_,length(y))
 if(nchar(s$bmi_var)) {e<-getd(x,s$bmi_file); j<-match(id,as.character(e[[grep('^bcsid$',names(e),ignore.case=TRUE,value=TRUE)[1]]]));bmi<-num(e[[s$bmi_var]][j])}
 add_summary('BCS70',s$age,stage(s$age),'adult_asthma_or_wheezy_bronchitis',y,bmi)
 q<-data.frame(id=id,y=y);names(q)[2]<-paste0('y',s$age);bad<-merge(bad,q,by='id',all.x=TRUE,sort=FALSE)
}
add_transition('BCS70','adult_asthma_or_wheezy_bronchitis',bad,c(38,42,46))
rm(x,a,b,bc,bad);gc()

cat('NextSteps\n');x<-readRDS(file.path(root,'NextSteps_full.rds'))
a<-getd(x,'ns8_2015_main_interview.dta');b<-getd(x,'ns9_2022_main_interview.dta')
y25<-yn(a$W8LOLR0Z,1,2);y32<-yn(b$W9LOLRAD,1,2)
add_summary('NextSteps',25,stage(25),'respiratory_impact_among_longstanding_illness',y25)
add_summary('NextSteps',32,stage(32),'respiratory_impact_among_longstanding_illness',y32)
add_summary('NextSteps',32,stage(32),'respiratory_condition_expected_12m',yn(b$W9DLOILCOND0A,1,2))
ns<-merge(data.frame(id=as.character(a$NSID),y25=y25),data.frame(id=as.character(b$NSID),y32=y32),by='id',all=TRUE,sort=FALSE)
add_transition('NextSteps','respiratory_impact_among_longstanding_illness',ns,c(25,32))
rm(x,a,b,ns);gc()

save_csv(do.call(rbind,summ),'cohort_wave_prevalence.csv')
save_csv(do.call(rbind,trans),'cohort_adjacent_transitions.csv')

# Restricted cubic spline is implemented as a natural cubic spline with fixed boundary knots.
# Age has four observed values, so use one internal knot; BMI has one internal knot.
cat('RCS\n')
d<-mcs_long[!is.na(mcs_long$wheeze12m)&is.finite(mcs_long$bmi),]
d<-d[order(d$id,d$age),];d$id<-factor(d$id)
save_csv(aggregate(wheeze12m~age,d,function(z)c(n=length(z),yes=sum(z))), 'MCS_rcs_sample_by_age_raw.csv')
library(splines)
library(geepack)
ak<-7; bk<-as.numeric(quantile(d$bmi,.5)); bounds<-as.numeric(quantile(d$bmi,c(.05,.95)))
d<-d[d$bmi>=bounds[1] & d$bmi<=bounds[2],]
af<-function(z) ns(z,knots=ak,Boundary.knots=c(3,11))
bf<-function(z) ns(z,knots=bk,Boundary.knots=bounds)
fit_add<-geeglm(wheeze12m~af(age)+bf(bmi),id=id,data=d,family=binomial,corstr='exchangeable')
fit_int<-geeglm(wheeze12m~af(age)*bf(bmi),id=id,data=d,family=binomial,corstr='exchangeable')
co<-coef(fit_int);vv<-vcov(fit_int); ix<-grep(':',names(co),fixed=TRUE)
W<-as.numeric(t(co[ix])%*%solve(vv[ix,ix,drop=FALSE],co[ix]));p<-pchisq(W,length(ix),lower.tail=FALSE)
save_csv(data.frame(term=names(co),estimate=as.numeric(co),se=sqrt(diag(vv)),z=as.numeric(co)/sqrt(diag(vv)),p=2*pnorm(-abs(as.numeric(co)/sqrt(diag(vv))))),'MCS_age_BMI_RCS_interaction_coefficients.csv')
save_csv(data.frame(model=c('additive','interaction'),n=nrow(d),children=length(unique(d$id)),age_knots='3|7|11',bmi_knots=paste(round(c(bounds[1],bk,bounds[2]),3),collapse='|'),wald_interaction_chisq=c(NA,W),interaction_df=c(NA,length(ix)),interaction_p=c(NA,p)), 'MCS_age_BMI_RCS_model_summary.csv')
pred_grid<-expand.grid(age=seq(3,11,length.out=81),bmi=as.numeric(quantile(d$bmi,c(.1,.5,.9))))
X<-model.matrix(delete.response(terms(fit_int)),pred_grid)
eta<-as.vector(X%*%co);se<-sqrt(rowSums((X%*%vv)*X))
pred_grid$prob<-plogis(eta);pred_grid$lo<-plogis(eta-1.96*se);pred_grid$hi<-plogis(eta+1.96*se)
save_csv(pred_grid,'MCS_age_BMI_RCS_predictions.csv')
saveRDS(list(additive=fit_add,interaction=fit_int,knots=list(age=c(3,ak,11),bmi=c(bounds[1],bk,bounds[2]))),file.path(out,'MCS_age_BMI_RCS_models.rds'),compress='xz')
cat('interaction W=',W,'p=',p,'n=',nrow(d),'children=',length(unique(d$id)),'\n')
