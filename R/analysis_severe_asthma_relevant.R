options(stringsAsFactors=FALSE, survey.lonely.psu='adjust')
suppressPackageStartupMessages({library(survey);library(sandwich);library(ggplot2)})
set.seed(20260928)
# ASCII staging directory avoids the local R/Windows CSV pathname encoding issue.
# Completed outputs are mirrored to the requested P1 folder by the PowerShell runner.
out <- 'E:/cohort_rds_stage_20260925/severe_asthma_relevant_analysis_20260928'
root <- 'E:/公共数据库/Four_Cohorts_Full_RDS_20260925'
v3 <- 'E:/cohort_rds_stage_20260925/life_course_v3_20260928'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
num <- function(x) suppressWarnings(as.numeric(unclass(x)))
yn <- function(x,yes=1,no=2){z<-num(x);ifelse(z==yes,1L,ifelse(z==no,0L,NA_integer_))}
getd <- function(b,key)b$datasets[[match(key,b$source_inventory$key)]]
write <- function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,fileEncoding='UTF-8',na='')
log_lines <- character()
qc <- function(s){log_lines<<-c(log_lines,s);cat(s,'\n')}
empty <- function(name,columns,note){z<-as.data.frame(setNames(replicate(length(columns),character(),simplify=FALSE),columns));write(z,name);qc(paste(name,note))}
ci_rows <- function(f,cohort,analysis,outcome,terms,type='RR',events=NA_integer_,n=NA_integer_){
  b<-coef(f);V<-vcov(f);terms<-intersect(terms,names(b));if(!length(terms))return(data.frame())
  s<-sqrt(diag(V)[terms]);v<-b[terms];data.frame(cohort,analysis,outcome,term=terms,n,events,scale=type,estimate=if(type=='RR')exp(v)else v,lower=if(type=='RR')exp(v-1.96*s)else v-1.96*s,upper=if(type=='RR')exp(v+1.96*s)else v+1.96*s,p=2*pnorm(-abs(v/s)),se=s)
}
contrast <- function(f,t1,t0,type='RR'){
  b<-coef(f);V<-vcov(f);stopifnot(t1%in%names(b),t0%in%names(b));v<-unname(b[t1]-b[t0]);s<-unname(sqrt(V[t1,t1]+V[t0,t0]-2*V[t1,t0]));c(estimate=if(type=='RR')exp(v)else v,lower=if(type=='RR')exp(v-1.96*s)else v-1.96*s,upper=if(type=='RR')exp(v+1.96*s)else v+1.96*s,p=2*pnorm(-abs(v/s)),se=s)
}
ess <- function(w)sum(w)^2/sum(w^2)
status <- list();flow <- list();diagnostics <- list();audit <- list();sens <- list();models <- list();family <- list()
add_flow<-function(cohort,analysis,stage,n,events=NA_integer_,note=''){flow[[length(flow)+1L]]<<-data.frame(cohort,analysis,stage,n,events,note)}
add_diag<-function(cohort,analysis,n,events,weight,minpred,maxpred,above1,notes='',parameters=NA_integer_){diagnostics[[length(diagnostics)+1L]]<<-data.frame(cohort,analysis,n,events,ESS=if(length(weight))ess(weight)else NA_real_,minimum_fitted=minpred,maximum_fitted=maxpred,n_fitted_over_one=above1,parameters=parameters,EPV=if(is.na(events)||is.na(parameters))NA_real_ else events/max(1,parameters),notes)}
add_audit<-function(cohort,wave,source_key,source_vars,derived,positive_code,negative_code,missing_rule,interpretation){audit[[length(audit)+1L]]<<-data.frame(cohort,wave,source_key,source_vars,derived,positive_code,negative_code,missing_rule,interpretation)}

# Figure contract: R quantitative grid; four figures show descriptive symptom severity,
# the lagged history gradient, the limited hospital data, and objective lung function.
# All estimates are observational associations. White background, editable SVG/PDF,
# PNG 400 dpi; source tables are the CSV files written here.

v <- readRDS(file.path(v3,'v3_cohort_datasets.rds'))

# ---- A. MCS severity derivation with one report selected per child and wave ----
m <- readRDS(file.path(root,'MCS_full.rds'))
make_wave <- function(prefix,key,age){
  d<-getd(m,key);id<-paste(as.character(d$MCSID),num(d[[paste0(prefix,'CNUM00')]]),sep=':')
  pnum<-num(d[[paste0(prefix,'PNUM00')]])
  ever<-yn(d[[paste0(prefix,'PASMA00')]]);ew<-yn(d[[paste0(prefix,'PWHEE00')]]);recent<-yn(d[[paste0(prefix,'PWHLY00')]])
  recent[is.na(recent)&ew==0]<-0L
  a<-num(d[[paste0(prefix,'PWHAN00')]])
  attack<-if(prefix=='B')ifelse(!is.na(a)&a>=0&a<99,as.integer(a>=4),NA_integer_) else if(prefix=='E')ifelse(a%in%1:4,as.integer(a%in%2:3),NA_integer_) else ifelse(a%in%1:4,as.integer(a%in%3:4),NA_integer_)
  sl<-num(d[[paste0(prefix,'PWHSD00')]]);sleep<-ifelse(sl%in%1:3,as.integer(sl==3),NA_integer_)
  speech<-yn(d[[paste0(prefix,'PWHSL00')]])
  z<-data.frame(id,age,parent_number=pnum,ever_asthma=ever,recent_wheeze=recent,attack4=attack,sleep_weekly=sleep,speech_limit=speech,ever_asthma_raw=num(d[[paste0(prefix,'PASMA00')]]),recent_wheeze_raw=num(d[[paste0(prefix,'PWHLY00')]]),attack4_raw=a,sleep_weekly_raw=sl,speech_limit_raw=num(d[[paste0(prefix,'PWHSL00')]]))
  z$completeness<-rowSums(!is.na(z[c('ever_asthma','recent_wheeze','attack4','sleep_weekly','speech_limit')]))
  duplicate_n<-sum(duplicated(z$id));qc(paste('MCS age',age,'raw rows',nrow(z),'duplicate child reports',duplicate_n))
  z<-z[order(z$id,-z$completeness,z$parent_number!=1),];z<-z[!duplicated(z$id),];stopifnot(!anyDuplicated(z$id))
  z$current_asthma_proxy<-ifelse(z$ever_asthma==1&z$recent_wheeze==1,1L,ifelse(z$ever_asthma==0|z$recent_wheeze==0,0L,NA_integer_))
  mat<-as.matrix(z[c('attack4','sleep_weekly','speech_limit')]);positive<-rowSums(mat==1,na.rm=TRUE)>0;allknown<-rowSums(!is.na(mat))==3
  severe<-ifelse(z$recent_wheeze==0,0L,ifelse(z$recent_wheeze==1&positive,1L,ifelse(z$recent_wheeze==1&allknown,0L,NA_integer_)))
  z$severe_wheeze_proxy<-severe
  z$asthma_background_severe<-ifelse(z$current_asthma_proxy==1&severe==1,1L,ifelse(z$current_asthma_proxy==0|severe==0,0L,NA_integer_))
  positives<-which(z$severe_wheeze_proxy==1)
  stopifnot(all(rowSums(mat[positives,,drop=FALSE]==1,na.rm=TRUE)>=1),all(z$recent_wheeze[positives]==1))
  add_audit('MCS',age,key,paste0(prefix,'PWHLY00; ',prefix,'PASMA00; ',prefix,'PWHAN00; ',prefix,'PWHSD00; ',prefix,'PWHSL00'),'recent/current-asthma/severe-wheeze','wheeze=yes; >=4 attacks OR weekly waking OR speech limitation','answered no wheeze is non-severe','99 and item skips remain missing; no recent wheeze=0 by definition','parent report; symptom proxy, not clinical severity')
  z
}
waves<-list(`3`=make_wave('B','MCS_00023',3),`5`=make_wave('C','MCS_00047',5),`7`=make_wave('D','MCS_00062',7),`11`=make_wave('E','MCS_00085',11))
rm(m);gc()
components<-c('recent_wheeze','ever_asthma','current_asthma_proxy','attack4','sleep_weekly','speech_limit','severe_wheeze_proxy','asthma_background_severe')
count_rows<-list();i<-0L
for(age in names(waves)){
 z<-waves[[age]];add_flow('MCS',paste0('age_',age),'unique_child_records',nrow(z));
 for(q in components){x<-z[[q]];raw<-z[[paste0(q,'_raw')]];i<-i+1L;count_rows[[i]]<-data.frame(age=as.integer(age),component=q,raw_unique_children=nrow(z),evaluable=sum(!is.na(x)),positive=sum(x==1,na.rm=TRUE),negative=sum(x==0,na.rm=TRUE),missing=sum(is.na(x)),routing_skip=if(length(raw))sum(raw==-1,na.rm=TRUE)else NA_integer_,other_raw_missing=if(length(raw))sum(is.na(x)&(is.na(raw)|raw!=-1),na.rm=TRUE)else NA_integer_,routing_note=if(q%in%c('attack4','sleep_weekly','speech_limit'))'component only asked on wheeze route; skips retained missing'else'negative recent wheeze is negative for composite')}
}
counts<-do.call(rbind,count_rows);write(counts,'02_MCS_severe_wheeze_counts.csv')

mc<-v$MCS
for(age in names(waves)){z<-waves[[age]][,c('id','recent_wheeze','current_asthma_proxy','attack4','sleep_weekly','speech_limit','severe_wheeze_proxy','asthma_background_severe')];names(z)[-1]<-paste0(names(z)[-1],'_',age);mc<-merge(mc,z,by='id',all.x=TRUE,sort=FALSE)}
stopifnot(!anyDuplicated(mc$id))
for(age in c(3,5,7,11)){oldy<-mc[[paste0('y',age)]];newy<-mc[[paste0('recent_wheeze_',age)]];stopifnot(sum(oldy!=newy,na.rm=TRUE)==0);qc(paste('MCS recent-wheeze agreement with V3 age',age,'discordant',sum(oldy!=newy,na.rm=TRUE)))}
qc(paste('MCS V3 baseline',nrow(mc),'linked age-11 severe outcome',sum(!is.na(mc$severe_wheeze_proxy_11)),'events',sum(mc$severe_wheeze_proxy_11==1,na.rm=TRUE)))
saveRDS(mc,file.path(out,'MCS_derived_analysis_frame.rds'),compress='xz')

covars<-c('male','low_birthweight','ses','pregnancy_continued_smoke','respondent_asthma','atopy3','bmi3_z','ever_asthma3')
get_weight<-function(d,age){w<-d[[paste0('weight_age',age)]];w[!is.finite(w)|w<=0]<-NA_real_;w}
fit_svy<-function(d,formula,weight){d$analysis_weight<-weight;design<-svydesign(ids=~family,strata=~design_stratum,weights=~analysis_weight,data=d,nest=TRUE);svyglm(formula,design=design,family=quasipoisson('log'))}
risk_results<-list();stage_results<-list()
for(age in c(5,7,11)){
 y<-paste0('severe_wheeze_proxy_',age);w<-get_weight(mc,age)
 z<-mc[complete.cases(mc[c(covars,'family','design_stratum',y)])&!is.na(w),];w<-w[match(z$id,mc$id)]
 add_flow('MCS',paste0('risk_factors_age_',age),'weighted_complete_case',nrow(z),sum(z[[y]]),paste('raw eligible',sum(!is.na(mc[[y]]))))
 if(sum(z[[y]])>=50){f<-fit_svy(z,reformulate(covars,y),w);risk_results[[as.character(age)]]<-ci_rows(f,'MCS',paste0('risk_factors_age_',age),'severe_wheeze_proxy',covars,n=nrow(z),events=sum(z[[y]]));mu<-fitted(f);add_diag('MCS',paste0('risk_factors_age_',age),nrow(z),sum(z[[y]]),w,min(mu),max(mu),sum(mu>1),'survey-weighted modified Poisson; family cluster and stratum',length(coef(f))-1)}
}
if(length(risk_results))risk_results<-do.call(rbind,risk_results) else risk_results<-data.frame()
write(risk_results,'03_MCS_severe_wheeze_models.csv')

# Age-11 category-specific risks: mild wheeze versus severe wheeze, each versus all children.
w<-get_weight(mc,11);y<-mc$severe_wheeze_proxy_11;rw<-mc$recent_wheeze_11
sel<-complete.cases(mc[c(covars,'family','design_stratum')])&!is.na(w)&!is.na(y)&!is.na(rw)
z<-mc[sel,];w<-w[sel];z$mild<-as.integer(z$recent_wheeze_11==1&z$severe_wheeze_proxy_11==0);z$severe<-as.integer(z$severe_wheeze_proxy_11==1)
stack<-rbind(transform(z,outcome=mild,severity='mild'),transform(z,outcome=severe,severity='severe'))
stack$severity<-factor(stack$severity,levels=c('mild','severe'));stack$analysis_weight<-rep(w,2)
des<-svydesign(ids=~family,strata=~design_stratum,weights=~analysis_weight,data=stack,nest=TRUE)
fs<-svyglm(reformulate(paste0('severity*(',paste(covars,collapse='+'),')'),'outcome'),design=des,family=quasipoisson('log'))
for(term in covars){it<-grep(paste0('^severitysevere:',term,'$|^',term,':severitysevere$'),names(coef(fs)),value=TRUE);if(length(it)==1){bb<-coef(fs)[it];se<-sqrt(vcov(fs)[it,it]);stage_results[[term]]<-data.frame(term,RR_ratio_severe_to_mild=exp(bb),lower=exp(bb-1.96*se),upper=exp(bb+1.96*se),p_interaction=2*pnorm(-abs(bb/se)),n=nrow(z),mild_events=sum(z$mild),severe_events=sum(z$severe))}}
stage_results<-do.call(rbind,stage_results);write(stage_results,'MCS_severity_architecture_interactions.csv')

# ---- B. Prior 3/5/7-year wheeze burden -> age-11 severity ----
mc$prior_count<-ifelse(complete.cases(mc[c('recent_wheeze_3','recent_wheeze_5','recent_wheeze_7')]),mc$recent_wheeze_3+mc$recent_wheeze_5+mc$recent_wheeze_7,NA_real_)
mc$prior_group<-factor(ifelse(is.na(mc$prior_count),NA_character_,ifelse(mc$prior_count==0,'never',ifelse(mc$prior_count==1,'isolated','repeated'))),levels=c('never','isolated','repeated'))
# Age-3 asthma is downstream of the start of this cumulative history and is omitted.
adj<-setdiff(covars,'ever_asthma3');ww<-get_weight(mc,11)
base<-mc[complete.cases(mc[c(adj,'family','design_stratum','severe_wheeze_proxy_11')])&!is.na(ww),];base$obs_history<-as.integer(!is.na(base$prior_group))
add_flow('MCS','history_to_age11','age11_outcome_and_covariates',nrow(base),sum(base$severe_wheeze_proxy_11==1),paste('history complete',sum(base$obs_history)))
if(length(unique(base$obs_history))==2){resp<-glm(reformulate(adj,'obs_history'),data=base,family=binomial());p<-pmin(.99,pmax(.05,predict(resp,type='response')));base$history_ipw<-mean(base$obs_history)/p}else base$history_ipw<-1
base$history_ipw<-pmin(quantile(base$history_ipw,.99),pmax(quantile(base$history_ipw,.01),base$history_ipw))
z<-base[base$obs_history==1,];w<-get_weight(z,11)*z$history_ipw
stopifnot(all(complete.cases(z[c('prior_group','severe_wheeze_proxy_11')])));qc(paste('MCS history sample',nrow(z),'severe age11',sum(z$severe_wheeze_proxy_11),'ESS',round(ess(w))))
add_flow('MCS','history_to_age11','history_complete_weighted',nrow(z),sum(z$severe_wheeze_proxy_11),paste('effective sample size',round(ess(w),1)))
hist_tab<-do.call(rbind,lapply(levels(z$prior_group),function(g){s<-z$prior_group==g;data.frame(cohort='MCS',age=11,prior_group=g,n=sum(s),severe_events=sum(z$severe_wheeze_proxy_11[s]),crude_risk=mean(z$severe_wheeze_proxy_11[s]),weighted_risk=weighted.mean(z$severe_wheeze_proxy_11[s],w[s]),ESS=ess(w[s]))}))
write(hist_tab,'MCS_prior_history_risks.csv')
f<-fit_svy(z,reformulate(c('prior_group',adj),'severe_wheeze_proxy_11'),w)
res<-ci_rows(f,'MCS','prior_wheeze_to_age11','severe_wheeze_proxy',c('prior_groupisolated','prior_grouprepeated'),n=nrow(z),events=sum(z$severe_wheeze_proxy_11))
history_ct<-contrast(f,'prior_grouprepeated','prior_groupisolated');res<-rbind(res,data.frame(cohort='MCS',analysis='prior_wheeze_to_age11',outcome='severe_wheeze_proxy',term='repeated_vs_isolated',n=nrow(z),events=sum(z$severe_wheeze_proxy_11),scale='RR',estimate=history_ct['estimate'],lower=history_ct['lower'],upper=history_ct['upper'],p=history_ct['p'],se=history_ct['se']))
models[['MCS_history']]<-res
add_diag('MCS','prior_wheeze_to_age11',nrow(z),sum(z$severe_wheeze_proxy_11),w,min(fitted(f)),max(fitted(f)),sum(fitted(f)>1),'age11 survey weight x history-completeness IPW, 1st/99th truncation',length(coef(f))-1)
# Empirical marginal standardization; model-conditional delta method.
ar<-list();V<-vcov(f);for(g in levels(z$prior_group)){nd<-z;nd$prior_group<-factor(g,levels=levels(z$prior_group));X<-model.matrix(formula(f),nd);mu<-as.numeric(exp(X%*%coef(f)));pw<-w/sum(w);risk<-sum(pw*mu);grad<-colSums(X*as.numeric(pw*mu));se<-sqrt(as.numeric(t(grad)%*%V%*%grad));ar[[g]]<-data.frame(prior_group=g,n=nrow(z),adjusted_risk=risk,lower=risk-1.96*se,upper=risk+1.96*se);attr(ar[[g]],'gradient')<-grad}
absrisk<-do.call(rbind,ar);write(absrisk,'MCS_prior_history_standardized_risks.csv')
rd<-list();for(g in c('isolated','never')){delta<-ar[['repeated']]$adjusted_risk-ar[[g]]$adjusted_risk;grad<-attr(ar[['repeated']],'gradient')-attr(ar[[g]],'gradient');se<-sqrt(as.numeric(t(grad)%*%V%*%grad));rd[[g]]<-data.frame(comparison=paste0('repeated_vs_',g),risk_difference=delta,lower=delta-1.96*se,upper=delta+1.96*se,p=2*pnorm(-abs(delta/se)))};write(do.call(rbind,rd),'MCS_standardized_risk_differences.csv')
# Weighted and unweighted, component-specific, asthma-background sensitivities.
for(method in c('weighted','unweighted'))for(outcome in c('severe_wheeze_proxy_11','attack4_11','sleep_weekly_11','speech_limit_11','asthma_background_severe_11')){
 q<-z[!is.na(z[[outcome]]),];wq<-if(method=='weighted')w[!is.na(z[[outcome]])]else rep(1,nrow(q));if(sum(q[[outcome]])<20)next
 fm<-fit_svy(q,reformulate(c('prior_group',adj),outcome),wq);rr<-ci_rows(fm,'MCS',paste0('history_',method),outcome,c('prior_groupisolated','prior_grouprepeated'),n=nrow(q),events=sum(q[[outcome]]));rr$method<-method;sens[[paste(method,outcome)]]<-rr
}
write(res,'04_childhood_burden_to_severe_wheeze.csv')
add_mcs_sens<-function(q,weight,outcome,method){
 if(nrow(q)<100||sum(q[[outcome]])<20)return(NULL)
 fm<-fit_svy(q,reformulate(c('prior_group',adj),outcome),weight);rr<-ci_rows(fm,'MCS','prior_history_sensitivity',outcome,c('prior_groupisolated','prior_grouprepeated'),n=nrow(q),events=sum(q[[outcome]]));dc<-contrast(fm,'prior_grouprepeated','prior_groupisolated');rr<-rbind(rr,data.frame(cohort='MCS',analysis='prior_history_sensitivity',outcome,term='repeated_vs_isolated',n=nrow(q),events=sum(q[[outcome]]),scale='RR',estimate=dc['estimate'],lower=dc['lower'],upper=dc['upper'],p=dc['p'],se=dc['se']));rr$method<-method;rr
}
sens[['MCS_survey_only']]<-add_mcs_sens(z,get_weight(z,11),'severe_wheeze_proxy_11','survey_weight_without_history_IPW')
sens[['MCS_unweighted_direct']]<-add_mcs_sens(z,rep(1,nrow(z)),'severe_wheeze_proxy_11','unweighted')
for(comp in c('attack4','sleep_weekly','speech_limit')){q<-z;ep<-paste0(comp,'_full_population');raw<-q[[paste0(comp,'_11')]];q[[ep]]<-ifelse(q$recent_wheeze_11==0,0L,raw);ok<-!is.na(q[[ep]]);sens[[ep]]<-add_mcs_sens(q[ok,],w[ok],ep,'survey_history_IPW_full_population')}
ok<-which(z$current_asthma_proxy_11==1);sens[['MCS_asthma_restricted']]<-add_mcs_sens(z[ok,],w[ok],'severe_wheeze_proxy_11','restricted_current_asthma_proxy')

# ---- C. NCDS age-33 conditional hospitalisation and BCS70 descriptive hospital data ----
ncds<-v$NCDS;ncds$child_group<-factor(ifelse(is.na(ncds$y11)|is.na(ncds$ch16),NA_character_,ifelse(ncds$y11+ncds$ch16==0,'never',ifelse(ncds$y11+ncds$ch16==1,'isolated','repeated'))),levels=c('never','isolated','repeated'))
n<-readRDS(file.path(root,'NCDS_full.rds'));a<-getd(n,'NCDS_00013');id<-as.character(a$ncdsid);stopifnot(!anyDuplicated(id));j<-match(ncds$id,id);ncds$admit33<-yn(a$n504025[j]);ncds$admit_count33<-num(a$n504026[j]);
qc(paste('NCDS n504025=1',sum(ncds$admit33==1,na.rm=TRUE),'n504025=0',sum(ncds$admit33==0,na.rm=TRUE),'count despite no',sum(ncds$admit33==0&!is.na(ncds$admit_count33),na.rm=TRUE),'count despite missing',sum(is.na(ncds$admit33)&!is.na(ncds$admit_count33))))
add_audit('NCDS',33,'NCDS_00013','n504025; n504026','admit33','1=yes','2=no','all other codes and non-routed records missing','past-12-month wheeze/asthma admission, conditional questionnaire route')
nh<-ncds[!is.na(ncds$admit33)&!is.na(ncds$child_group),];htab<-do.call(rbind,lapply(levels(ncds$child_group),function(g){q<-nh[nh$child_group==g,];data.frame(cohort='NCDS',endpoint='age33_wheeze_asthma_admission',exposure=g,n=nrow(q),events=sum(q$admit33),crude_risk=if(nrow(q))mean(q$admit33)else NA_real_,comparison='within answered/eligible route',model_status='descriptive_only')}))
add_flow('NCDS','age33_admission','full_cohort',nrow(ncds));add_flow('NCDS','age33_admission','outcome_answered',sum(!is.na(ncds$admit33)),sum(ncds$admit33==1,na.rm=TRUE));add_flow('NCDS','age33_admission','child_history_and_outcome',nrow(nh),sum(nh$admit33))
min_ev<-min(htab$events);qc(paste('NCDS age33 analysed admissions',sum(nh$admit33),'group minimum',min_ev,'=> no multivariable RR due to sparse event cells and selective route'))
rm(n,a);gc()
b<-readRDS(file.path(root,'BCS70_full.rds'));bcs<-v$BCS70
da<-getd(b,'BCS70_00011');db<-getd(b,'BCS70_00006');id10<-as.character(da$bcsid);id16<-as.character(db$bcsid);stopifnot(!anyDuplicated(id10),!anyDuplicated(id16));j10<-match(bcs$id,id10);j16<-match(bcs$id,id16)
bcs$admit10_positive<-ifelse(num(da$b7_18[j10])==1,1L,NA_integer_)
bcs$admit16_attack_positive<-ifelse(num(db$ob9_14[j16])==1,1L,NA_integer_)
bcs$admit16_since10_positive<-ifelse(num(db$ob13_25[j16])%in%c(1,3),1L,NA_integer_)
bcs$child_group<-factor(ifelse(is.na(bcs$y10)|is.na(bcs$ch16),NA_character_,ifelse(bcs$y10+bcs$ch16==0,'never',ifelse(bcs$y10+bcs$ch16==1,'isolated','repeated'))),levels=c('never','isolated','repeated'))
bh<-do.call(rbind,lapply(c('admit10_positive','admit16_attack_positive','admit16_since10_positive'),function(ep)do.call(rbind,lapply(levels(bcs$child_group),function(g){q<-bcs[bcs$child_group==g,];data.frame(cohort='BCS70',endpoint=ep,exposure=g,n=sum(!is.na(q[[ep]])),events=sum(q[[ep]]==1,na.rm=TRUE),crude_risk=NA_real_,comparison='positive records only; no valid negative denominator',model_status='descriptive_only')}))))
add_audit('BCS70',10,'BCS70_00011','b7_18','admit10_positive','1=yes','no verified negative code','blank/skip is unknown','contemporaneous positive-only wheeze admission marker')
add_audit('BCS70',16,'BCS70_00006','ob9_14; ob13_25','admit16_positive','ob9_14=1; ob13_25 in 1,3','none verified','2 and 4 are no response; skips unknown','historical windows; descriptive only')
add_flow('BCS70','hospital','age10_positive',sum(bcs$admit10_positive==1,na.rm=TRUE));add_flow('BCS70','hospital','age16_attack_positive',sum(bcs$admit16_attack_positive==1,na.rm=TRUE));add_flow('BCS70','hospital','age16_since10_positive',sum(bcs$admit16_since10_positive==1,na.rm=TRUE))
hospital<-rbind(htab,bh);write(hospital,'05_hospitalisation_analysis.csv');rm(b,da,db);gc()

# ---- D. NCDS objective midlife lung function ----
n<-readRDS(file.path(root,'NCDS_full.rds'));lf<-getd(n,'NCDS_00046');id<-as.character(lf$ncdsid);stopifnot(!anyDuplicated(id));j<-match(ncds$id,id)
ncds$age_lf<-num(lf$page[j]);ncds$height_cm<-num(lf$htres[j]);ncds$FEV1_L<-num(lf$htfev[j]);ncds$FVC_L<-num(lf$htfvc[j]);ncds$PEF_Lmin<-num(lf$htpf[j]);ncds$smoking_midlife<-num(lf$smok42_1[j]);ncds$lf_complete<-num(lf$lfoc[j]);ncds$three_blow<-num(lf$done3[j]);ncds$no_read_flag<-num(lf$noread[j]);ncds$inhaler_24h<-num(lf$binhaler[j])
gap<-function(d,stem){a<-sapply(c(stem,paste0(stem,2:5)),function(v)num(d[[v]]));apply(a,1,function(x){x<-sort(x[is.finite(x)&x>.2&x<9],decreasing=TRUE);if(length(x)>=2)x[1]-x[2]else NA_real_})}
ncds$FEV_repeat_gap<-gap(lf,'fev')[j];ncds$FVC_repeat_gap<-gap(lf,'fvc')[j]
ncds$weight_kg<-num(lf$bweight[j]);ncds$bmi_midlife<-ncds$weight_kg/(ncds$height_cm/100)^2;ncds$bmi_midlife[ncds$weight_kg<30|ncds$weight_kg>250|ncds$bmi_midlife<12|ncds$bmi_midlife>60]<-NA_real_
ncds$pair_screen<-!is.na(ncds$FEV1_L)&!is.na(ncds$FVC_L)&ncds$FEV1_L>.2&ncds$FEV1_L<8&ncds$FVC_L>.3&ncds$FVC_L<9&ncds$FEV1_L<=ncds$FVC_L&ncds$lf_complete==1&ncds$no_read_flag==2
ncds$strict_qc<-ncds$pair_screen&ncds$three_blow==1&is.finite(ncds$height_cm)&ncds$height_cm>=120&ncds$height_cm<=210&ncds$age_lf%in%44:46
ncds$pair_screen[is.na(ncds$pair_screen)]<-FALSE;ncds$strict_qc[is.na(ncds$strict_qc)]<-FALSE
ncds$ratio<-ifelse(ncds$pair_screen,ncds$FEV1_L/ncds$FVC_L,NA_real_)
ncds$smoking_midlife[!ncds$smoking_midlife%in%0:2]<-NA_real_
qc(paste('NCDS spirometry paired',sum(ncds$pair_screen,na.rm=TRUE),'strict QC three acceptable blows+height',sum(ncds$strict_qc,na.rm=TRUE)))
qc_tab<-data.frame(stage=c('biomedical_records','completed_lung_test','plausible_paired_FEV1_FVC','no_read_flag_clear','three_acceptable_blows','plausible_measured_height','strict_QC','strict_QC_with_child_burden','strict_QC_model_complete'),n=c(length(id),sum(num(lf$lfoc)==1,na.rm=TRUE),sum(ncds$pair_screen,na.rm=TRUE),sum(num(lf$noread)==2,na.rm=TRUE),sum(num(lf$done3)==1,na.rm=TRUE),sum(num(lf$htres)>=120&num(lf$htres)<=210,na.rm=TRUE),sum(ncds$strict_qc,na.rm=TRUE),sum(ncds$strict_qc&!is.na(ncds$child_group),na.rm=TRUE),sum(ncds$strict_qc&complete.cases(ncds[c('child_group','male','ses','low_birthweight','pregnancy_continued_smoke','smoking_midlife')]),na.rm=TRUE)),note=c('','at least one technically satisfactory blow','limits FEV1 .2-8 L, FVC .3-9 L, FEV1<=FVC','noread=2','done3=1','height 120-210 cm','primary QC','two recent childhood windows','includes model covariates'))
write(qc_tab,'06_NCDS_lung_function_QC.csv')
repeat_qc<-data.frame(metric=c('FEV_top_two_gap_le_150ml','FVC_top_two_gap_le_150ml','both_top_two_gap_le_150ml','inhaler_past24h_yes','inhaler_past24h_no'),n=c(sum(ncds$strict_qc&ncds$FEV_repeat_gap<=.15,na.rm=TRUE),sum(ncds$strict_qc&ncds$FVC_repeat_gap<=.15,na.rm=TRUE),sum(ncds$strict_qc&ncds$FEV_repeat_gap<=.15&ncds$FVC_repeat_gap<=.15,na.rm=TRUE),sum(ncds$inhaler_24h==1,na.rm=TRUE),sum(ncds$inhaler_24h==2,na.rm=TRUE)),note='raw blow repeatability screen; not a complete ATS quality grade; inhaler type unknown')
write(repeat_qc,'NCDS_repeated_blows_QC.csv')
add_audit('NCDS','44-46','NCDS_00046','htfev; htfvc; htpf; lfoc; noread; done3; htres; page; smok42_1','FEV1/FVC continuous','strict QC requires >=3 acceptable blows and valid pair','n/a','negative/sentinel values excluded','physiological indices, no bronchodilator reversibility')
add_flow('NCDS','lung_function','strict_QC',sum(ncds$strict_qc,na.rm=TRUE));
rm(n,lf);gc()
adjlf<-c('child_group','male','height_cm','age_lf','smoking_midlife','ses','low_birthweight','pregnancy_continued_smoke')
lz<-ncds[ncds$strict_qc&complete.cases(ncds[adjlf]),];lz$smoking_midlife<-factor(lz$smoking_midlife,levels=0:2);lz$age_lf<-as.numeric(lz$age_lf)
add_flow('NCDS','lung_function','model_complete',nrow(lz),NA_integer_,paste('child groups',paste(table(lz$child_group),collapse='/')))
lf_results<-list();lf_summary<-list();
for(outcome in c('FEV1_L','FVC_L','ratio')){
 f<-lm(reformulate(adjlf,outcome),data=lz);V<-vcovHC(f,type='HC3');b<-coef(f);for(term in c('child_groupisolated','child_grouprepeated')){se<-sqrt(V[term,term]);lf_results[[paste(outcome,term)]]<-data.frame(cohort='NCDS',endpoint=outcome,analysis='strict_QC_adjusted_linear',term,n=nrow(lz),mean=mean(lz[[outcome]]),estimate=b[term],lower=b[term]-1.96*se,upper=b[term]+1.96*se,p=2*pnorm(-abs(b[term]/se)),standardized_beta=b[term]/sd(lz[[outcome]]))}
 # Explicit robust contrast, since lm's ordinary vcov is not the HC3 covariance.
 dd<-b['child_grouprepeated']-b['child_groupisolated'];se<-sqrt(V['child_grouprepeated','child_grouprepeated']+V['child_groupisolated','child_groupisolated']-2*V['child_grouprepeated','child_groupisolated'])
 lf_results[[paste(outcome,'direct')]]<-data.frame(cohort='NCDS',endpoint=outcome,analysis='strict_QC_adjusted_linear',term='repeated_vs_isolated',n=nrow(lz),mean=mean(lz[[outcome]]),estimate=dd,lower=dd-1.96*se,upper=dd+1.96*se,p=2*pnorm(-abs(dd/se)),standardized_beta=dd/sd(lz[[outcome]]))
 lf_summary[[outcome]]<-do.call(rbind,lapply(levels(lz$child_group),function(g){q<-lz[lz$child_group==g,];data.frame(endpoint=outcome,child_group=g,n=nrow(q),mean=mean(q[[outcome]]),sd=sd(q[[outcome]]))}))
 add_diag('NCDS',paste0('lung_',outcome),nrow(lz),NA_integer_,numeric(),min(fitted(f)),max(fitted(f)),NA_integer_,paste('HC3 robust SE; R2',round(summary(f)$r.squared,3),'max Cook D',round(max(cooks.distance(f)),3)))
}
lf_results<-do.call(rbind,lf_results);write(lf_results,'07_NCDS_lung_function_models.csv');write(do.call(rbind,lf_summary),'NCDS_lung_function_by_child_group.csv')
cut<-quantile(lz$ratio,.1);lz$low_ratio<-as.integer(lz$ratio<=cut);f<-glm(low_ratio~child_group+male+height_cm+age_lf+smoking_midlife+ses+low_birthweight+pregnancy_continued_smoke,data=lz,family=poisson('log'));V<-vcovHC(f,type='HC3');for(term in c('child_groupisolated','child_grouprepeated')){bb<-coef(f)[term];se<-sqrt(V[term,term]);sens[[paste('lower_tail',term)]]<-data.frame(cohort='NCDS',analysis='lung_lower_tail',outcome='FEV1_FVC_bottom10pct',term,n=nrow(lz),events=sum(lz$low_ratio),scale='RR',estimate=exp(bb),lower=exp(bb-1.96*se),upper=exp(bb+1.96*se),p=2*pnorm(-abs(bb/se)),se=se,method='strict_QC')}
qc(paste('NCDS FEV1/FVC bottom-decile cut',round(cut,4),'events',sum(lz$low_ratio)))

# Response-IPW estimates are sensitivity analyses for MAR given childhood variables.
baseline_vars<-c('child_group','male','ses','low_birthweight','pregnancy_continued_smoke')
rb<-ncds[complete.cases(ncds[baseline_vars])&!ncds$dead42,];rb$observed_lf_model<-as.integer(rb$strict_qc&complete.cases(rb[adjlf]))
fr<-glm(reformulate(baseline_vars,'observed_lf_model'),data=rb,family=binomial());pr<-pmin(.99,pmax(.05,predict(fr,type='response')));rb$response_ipw<-mean(rb$observed_lf_model)/pr
qq<-quantile(rb$response_ipw,c(.01,.99));rb$response_ipw<-pmin(qq[2],pmax(qq[1],rb$response_ipw));lz$response_ipw<-rb$response_ipw[match(lz$id,rb$id)]
qc(paste('NCDS response-IPW baseline',nrow(rb),'observed',sum(rb$observed_lf_model),'ESS',round(ess(lz$response_ipw))))
linear_sensitivity<-function(q,outcome,method,weight=rep(1,nrow(q)),extra=character()){
 ad<-c(adjlf,extra);ok<-complete.cases(q[ad])&is.finite(q[[outcome]])&is.finite(weight)&weight>0;q<-q[ok,];weight<-weight[ok];if(nrow(q)<100)return(NULL)
 q$smoking_midlife<-factor(q$smoking_midlife,levels=0:2);q$.w<-weight
 fm<-lm(reformulate(ad,outcome),data=q,weights=.w);V<-vcovHC(fm,type='HC3');b<-coef(fm);r<-list()
 for(term in c('child_groupisolated','child_grouprepeated')){se<-sqrt(V[term,term]);r[[term]]<-data.frame(cohort='NCDS',analysis='lung_sensitivity',outcome,term,n=nrow(q),events=NA_integer_,scale='beta',estimate=unname(b[term]),lower=unname(b[term]-1.96*se),upper=unname(b[term]+1.96*se),p=2*pnorm(-abs(b[term]/se)),se,method)}
 dd<-unname(b['child_grouprepeated']-b['child_groupisolated']);se<-sqrt(V['child_grouprepeated','child_grouprepeated']+V['child_groupisolated','child_groupisolated']-2*V['child_grouprepeated','child_groupisolated']);r[['direct']]<-data.frame(cohort='NCDS',analysis='lung_sensitivity',outcome,term='repeated_vs_isolated',n=nrow(q),events=NA_integer_,scale='beta',estimate=dd,lower=dd-1.96*se,upper=dd+1.96*se,p=2*pnorm(-abs(dd/se)),se,method);do.call(rbind,r)
}
for(ep in c('FEV1_L','ratio')){
 sens[[paste0(ep,'_responseIPW')]]<-linear_sensitivity(lz,ep,'response_IPW',lz$response_ipw)
 sens[[paste0(ep,'_BMI')]]<-linear_sensitivity(lz,ep,'additional_midlife_BMI',extra='bmi_midlife')
 relaxed<-ncds[ncds$pair_screen&ncds$height_cm>=120&ncds$height_cm<=210&ncds$age_lf%in%44:46,];relaxed<-relaxed[!is.na(relaxed$id),];sens[[paste0(ep,'_relaxed')]]<-linear_sensitivity(relaxed,ep,'one_acceptable_blow_QC')
 q<-lz[which(lz$FEV_repeat_gap<=.15&lz$FVC_repeat_gap<=.15),];sens[[paste0(ep,'_repeatability')]]<-linear_sensitivity(q,ep,'raw_top2_gaps_le_150ml')
 q<-lz[which(lz$inhaler_24h==2),];sens[[paste0(ep,'_noinhaler')]]<-linear_sensitivity(q,ep,'no_inhaler_past24h')
 simple<-lm(reformulate(c('male','height_cm','age_lf'),ep),data=lz);sdr<-sd(residuals(simple));q<-lz;q$residual_z<-residuals(simple)/sdr
 rr<-linear_sensitivity(q,'residual_z','cohort_internal_age_sex_height_residual_z');rr$outcome<-paste0(ep,'_residual_z');sens[[paste0(ep,'_residual_z')]]<-rr
}
saveRDS(ncds,file.path(out,'NCDS_derived_analysis_frame.rds'),compress='xz')
saveRDS(lz,file.path(out,'NCDS_lung_model_frame.rds'),compress='xz')

# ---- Multiplicity, translation matrix, and quality artifacts ----
all_sens<-if(length(sens))do.call(rbind,lapply(sens,function(x){z<-x;for(q in setdiff(c('cohort','analysis','outcome','term','n','events','scale','estimate','lower','upper','p','se','method'),names(z)))z[[q]]<-NA;z[c('cohort','analysis','outcome','term','n','events','scale','estimate','lower','upper','p','se','method')]}))else data.frame()
write(all_sens,'08_sensitivity_results.csv')
families<-list(
 A_MCS_risk_factors=data.frame(family='A_MCS_risk_factors',analysis=risk_results$analysis,term=risk_results$term,p=risk_results$p),
 A_MCS_severity_architecture=data.frame(family='A_MCS_severity_architecture',analysis='age11_mild_vs_severe',term=stage_results$term,p=stage_results$p_interaction),
 B_child_history=data.frame(family='B_child_history',analysis=res$analysis,term=res$term,p=res$p),
 D_lung_function=data.frame(family='D_lung_function',analysis=lf_results$endpoint,term=lf_results$term,p=lf_results$p)
)
mult<-do.call(rbind,lapply(families,function(x){x$q_BH<-p.adjust(x$p,'BH');x}));write(mult,'09_multiple_testing.csv')
if(nrow(all_sens)){ss<-data.frame(family=ifelse(all_sens$cohort=='MCS','B_component_and_weight_sensitivity','D_lung_sensitivity'),analysis=paste(all_sens$method,all_sens$outcome,sep=':'),term=all_sens$term,p=all_sens$p);ss$q_BH<-ave(ss$p,ss$family,FUN=function(x)p.adjust(x,'BH'));mult<-rbind(mult,ss)}
mult<-rbind(mult,data.frame(family='C_hospitalisation',analysis='sparse_or_positive_only',term='not_estimable',p=NA_real_,q_BH=NA_real_));write(mult,'09_multiple_testing.csv')
write(do.call(rbind,flow),'10_sample_flow.csv');write(do.call(rbind,diagnostics),'11_model_diagnostics.csv');write(do.call(rbind,audit),'01_variable_derivation_audit.csv')
mat<-data.frame(
 cohort=c('MCS','NCDS','BCS70','NCDS','NCDS'),
 exposure_phenotype=c('prior age3/5/7 wheeze 0/1/>=2','age11/16 wheeze 0/1/2','age10/16 wheeze 0/1/2','age11/16 wheeze 0/1/2','age11/16 wheeze 0/1/2'),
 endpoint=c('severe-wheeze symptom proxy','age33 wheeze/asthma admission','historical wheeze/asthma admission','FEV1 (L)','FEV1/FVC'),age=c('11','33','10/16','44-46','44-46'),
 n=c(nrow(z),nrow(nh),NA_integer_,nrow(lz),nrow(lz)),
 event_n_or_mean=c(as.character(sum(z$severe_wheeze_proxy_11)),as.character(sum(nh$admit33)),paste('positive-only',sum(bcs$admit16_since10_positive==1,na.rm=TRUE)),as.character(round(mean(lz$FEV1_L),3)),as.character(round(mean(lz$ratio),3))),
 effect_estimate=c(history_ct['estimate'],NA_real_,NA_real_,lf_results$estimate[lf_results$endpoint=='FEV1_L'&lf_results$term=='child_grouprepeated'],lf_results$estimate[lf_results$endpoint=='ratio'&lf_results$term=='child_grouprepeated']),
 lower=c(history_ct['lower'],NA_real_,NA_real_,lf_results$lower[lf_results$endpoint=='FEV1_L'&lf_results$term=='child_grouprepeated'],lf_results$lower[lf_results$endpoint=='ratio'&lf_results$term=='child_grouprepeated']),
 upper=c(history_ct['upper'],NA_real_,NA_real_,lf_results$upper[lf_results$endpoint=='FEV1_L'&lf_results$term=='child_grouprepeated'],lf_results$upper[lf_results$endpoint=='ratio'&lf_results$term=='child_grouprepeated']),
 p=c(history_ct['p'],NA_real_,NA_real_,lf_results$p[lf_results$endpoint=='FEV1_L'&lf_results$term=='child_grouprepeated'],lf_results$p[lf_results$endpoint=='ratio'&lf_results$term=='child_grouprepeated']),
 interpretation=c('lagged symptom burden association','conditional questionnaire route; sparse events','positive-only codes; no negative risk denominator','objective physiology, adjusted association','objective physiology, adjusted association'),
 comparability_level=c('cohort-specific','cohort-specific','descriptive only','NCDS-only objective','NCDS-only objective'),
 placement=c('main-candidate','supplement','supplement','main-candidate','main-candidate'))
mat$comparison<-c('repeated versus isolated','descriptive counts by childhood group','descriptive positive records','repeated versus never','repeated versus never')
write(mat,'12_translation_evidence_matrix.csv')
risk_results$q_BH<-p.adjust(risk_results$p,'BH');stage_results$q_BH<-p.adjust(stage_results$p_interaction,'BH');res$q_BH<-p.adjust(res$p,'BH');lf_results$q_BH<-p.adjust(lf_results$p,'BH')
all_sens$q_BH<-ave(all_sens$p,all_sens$cohort,FUN=function(x)p.adjust(x,'BH'))
write(risk_results,'03_MCS_severe_wheeze_models.csv');write(stage_results,'MCS_severity_architecture_interactions.csv');write(res,'04_childhood_burden_to_severe_wheeze.csv');write(lf_results,'07_NCDS_lung_function_models.csv');write(all_sens,'08_sensitivity_results.csv')
saveRDS(mc,file.path(out,'MCS_derived_analysis_frame.rds'),compress='xz')
writeLines(log_lines,file.path(out,'quality_checks.txt'),useBytes=TRUE);writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'),useBytes=TRUE)
qc('Analysis tables complete; plotting stage follows')
