root<-'E:/公共数据库/Four_Cohorts_Full_RDS_20260925';out<-'E:/cohort_rds_stage_20260925/life_course_v3_20260928';dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw<-list();sink(file.path(out,'child_source_audit.txt'))
for(c in c('NCDS','BCS70')){
 x<-readRDS(file.path(root,paste0(c,'_full.rds')));fns<-basename(gsub('\\\\','/',x$source_inventory$source_entry));targets<-if(c=='NCDS')'ncds0123.dta'else c('sn3723.dta','bcs7016x.dta')
 for(fn in targets){i<-which(fns==fn);stopifnot(length(i)==1);d<-x$datasets[[i]];labs<-vapply(d,function(v){z<-attr(v,'label');if(is.null(z))''else as.character(z)[1]},'');jj<-grep('asthma|wheez|bronch|recent attack|ob9_|^n26[0-9] |^ncdsid |^bcsid ',paste(names(d),labs),ignore.case=TRUE);raw[[paste(c,fn,sep=':')]]<-d[unique(c(grep('^ncdsid$|^bcsid$',names(d),ignore.case=TRUE),jj))]
 cat('\n',c,fn,nrow(d),'\n');for(j in jj){v<-d[[j]];vl<-attr(v,'labels');tb<-sort(table(as.character(v),useNA='always'),decreasing=TRUE);cat(names(d)[j],labs[j],paste(names(vl),vl,sep=':',collapse='|'),'FREQ',paste(names(head(tb,15)),head(tb,15),sep=':',collapse='|'),'\n')}
 }
 rm(x);gc()
}
sink();saveRDS(raw,file.path(out,'child_selected_sources.rds'));cat('Audit complete\n')
