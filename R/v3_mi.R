library(mice)
t<-readLines('v3_models.R');ii<-grep('^if\\(Sys.getenv',t);eval(parse(text=t[seq_len(ii[1]-1)]))
set.seed(20260929)
pool<-function(fs,c,endpoint,terms){m<-length(fs);B<-sapply(fs,function(ff)coef(ff$f)[terms]);U<-sapply(fs,function(ff)diag(ff$V)[terms]);q<-rowMeans(B);ub<-rowMeans(U);between<-apply(B,1,var);total<-ub+(1+1/m)*between;df<-ifelse(between<1e-20,Inf,(m-1)*(1+ub/((1+1/m)*between))^2);crit<-qt(.975,df);data.frame(cohort=c,endpoint,method='MI_IPCW',term=terms,log_RR=q,SE=sqrt(total),RR=exp(q),lower=exp(q-crit*sqrt(total)),upper=exp(q+crit*sqrt(total)),p=2*pt(-abs(q/sqrt(total)),df),imputations=m,df,lambda=(1+1/m)*between/total,FMI=((1+1/m)*between/total)+(1-(1+1/m)*between/total)*2/(df+3),MCSE=sqrt(between/m),MCSE_over_SE=sqrt(between/m)/sqrt(total))}
rhat<-function(imp,vars){n<-dim(imp$chainMean)[2];its<-seq.int(max(1,n-9),n);sapply(vars,function(v){x<-imp$chainMean[v,its,,drop=FALSE];x<-matrix(x,nrow=length(its));W<-mean(apply(x,2,var));B<-length(its)*var(colMeans(x));if(!is.finite(W)||W==0)NA_real_ else sqrt(((length(its)-1)/length(its)*W+B/length(its))/W)})}
allresults<-list();mf<-list();conv<-list();traces<-list();individual<-list()
for(c in c('MCS','NCDS','BCS70')){
 cat('\nSTART MI',c,format(Sys.time()),'\n');flush.console()
 d<-data[[c]];vars<-if(c=='MCS')mcsvars else adultvars;ages<-switch(c,MCS=c(3,5,7,11),NCDS=c(23,33,42,50,55,62),BCS70=c(26,30,34,38,42,46));yv<-paste0('y',ages)
 if(c=='MCS')d<-d[!is.na(d$y3)&complete.cases(d[c('basewt','family','design_stratum')]),]else d<-d[rowSums(!is.na(d[c(vars,yv)]))>0,]
 impdata<-d[vars];for(v in setdiff(vars,'bmi3_z'))impdata[[v]]<-factor(impdata[[v]],levels=0:1)
 for(age in ages){y<-d[[paste0('y',age)]];impdata[[paste0('value',age)]]<-ifelse(is.na(y),0,y);obs<-as.integer(!is.na(y));if(length(unique(obs))>1)impdata[[paste0('observed',age)]]<-obs}
 if(c=='MCS'){impdata$log_basewt<-log(d$basewt);impdata$stratum<-factor(d$design_stratum)}
 meth<-setNames(rep('',ncol(impdata)),names(impdata));meth[vars]<-ifelse(vars=='bmi3_z','pmm','logreg');meth[vars[colSums(is.na(impdata[vars]))==0]]<-''
 pm<-make.predictorMatrix(impdata);diag(pm)<-0;pm[setdiff(names(impdata),vars),]<-0
 imp<-mice(impdata,m=20,maxit=15,method=meth,predictorMatrix=pm,seed=20261000+match(c,c('MCS','NCDS','BCS70')),printFlag=TRUE)
 rh<-rhat(imp,vars);if(any(rh>1.1,na.rm=TRUE)){cat('Extending chains by 15 iterations:',c,'\n');imp<-mice.mids(imp,maxit=15,printFlag=TRUE);rh<-rhat(imp,vars)}
 saveRDS(list(ids=d$id,mids=imp),file.path(out,paste0(c,'_MI_object.rds')))
 mf[[c]]<-data.frame(cohort=c,term=vars,eligible=nrow(d),missing=colSums(is.na(impdata[vars])),fraction=colMeans(is.na(impdata[vars])))
 conv[[c]]<-data.frame(cohort=c,term=vars,iterations=imp$iteration,Rhat_mean_last10=rh)
 if(!is.null(imp$loggedEvents))write.csv(imp$loggedEvents,file.path(out,paste0(c,'_MI_loggedEvents.csv')),row.names=FALSE)
 cm<-imp$chainMean;for(v in vars)for(j in 1:20)traces[[paste(c,v,j)]]<-data.frame(cohort=c,term=v,imputation=j,iteration=seq_len(dim(cm)[2]),chain_mean=cm[v,,j])
 fs0<-fs1<-fsa<-list();flows<-list()
 for(j in 1:20){cat('Fitting',c,j,'\n');flush.console();completed<-complete(imp,j);u<-d;for(v in vars)u[[v]]<-if(v=='bmi3_z')completed[[v]]else as.numeric(as.character(completed[[v]]));stopifnot(!anyNA(u[vars]));ll<-ipcw_long(u,ages,vars,c,c=='MCS');z<-ll$monotone
  if(c=='MCS'){z<-droplevels(z[z$age>3,]);for(origin in 0:1){q<-z[z$previous==origin,];ff<-estimate(q,reformulate(c('wave',vars),'y'),c);if(origin==0)fs0[[j]]<-ff else fs1[[j]]<-ff;individual[[paste(c,origin,j)]]<-cbind(imputation=j,extract(ff,c,if(origin==0)'transition_0_to_1'else'transition_1_to_1','MI_IPCW',vars))}}else{ff<-estimate(z,reformulate(c('wave',vars),'y'),c);fsa[[j]]<-ff;individual[[paste(c,j)]]<-cbind(imputation=j,extract(ff,c,'adult_repeated_burden','MI_IPCW',vars))}
  flows[[j]]<-data.frame(cohort=c,imputation=j,eligible=nrow(d),people=length(unique(z$id)),rows=nrow(z),events=sum(z$y),final_people=sum(z$age==max(ages)),max_weight=max(z$w))
 }
 if(c=='MCS'){allresults[[paste(c,0)]]<-pool(fs0,c,'transition_0_to_1',vars);allresults[[paste(c,1)]]<-pool(fs1,c,'transition_1_to_1',vars)}else allresults[[c]]<-pool(fsa,c,'adult_repeated_burden',vars)
 write.csv(do.call(rbind,flows),file.path(out,paste0(c,'_MI_sample_flow.csv')),row.names=FALSE)
 write.csv(do.call(rbind,allresults),file.path(out,'MI_IPCW_core_results.csv'),row.names=FALSE)
 write.csv(do.call(rbind,mf),file.path(out,'MI_missingness.csv'),row.names=FALSE);write.csv(do.call(rbind,conv),file.path(out,'MI_chain_diagnostics.csv'),row.names=FALSE);write.csv(do.call(rbind,traces),file.path(out,'MI_chain_traces.csv'),row.names=FALSE)
 cat('DONE MI',c,format(Sys.time()),'\n');flush.console()
}
write.csv(do.call(rbind,individual),file.path(out,'MI_individual_fit_results.csv'),row.names=FALSE);cat('MI+IPCW all core models complete\n')
