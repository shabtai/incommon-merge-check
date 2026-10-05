suppressMessages(library(survival))
pub<-readRDS('/private/tmp/incommon-survival.rds'); pub<-pub[pub$fit_type=='general',]
src<-read.delim('/private/tmp/incommon-paper-sample.txt',comment.char='#',check.names=FALSE)
keys<-list(c('PAAD','KRAS','PAAD'),c('PAAD','TP53','PAAD'),c('OV','TP53','HGSOC'),c('MEL','NRAS','SKCM'),c('MEL','TP53','SKCM'),c('HCC','TP53','HCC'),c('CRC:CRC MSS','BRAF:p.V600E','COAD'))
set.seed(1); out<-list()
for(k in keys){
 i<-which(pub$tumor_type==k[1]&pub$gene==k[2]); p<-pub$cox_fit[[i]]; e<-environment(p$terms); d<-eval(p$call$data,envir=e)
 d$orig<-src$ONCOTREE_CODE[match(d$sample,src$SAMPLE_ID)]
 used<-complete.cases(model.frame(formula(p),data=d,na.action=na.pass)); d<-d[used,]
 keep<-sum(d$orig==k[3]); terms<-grep('^class',names(coef(p)),value=TRUE)
 sims<-t(replicate(300,{s<-d[sample(nrow(d),keep),]; f<-try(coxph(formula(p),data=s,ties=p$method),silent=TRUE)
   if(inherits(f,'try-error')) rep(NA,2*length(terms)) else {sm<-summary(f)$coef; c(exp(sm[terms,1]),sm[terms,5])}}))
 nt<-length(terms)
 for(j in seq_len(nt)) out[[length(out)+1]]<-data.frame(tt=k[1],gene=k[2],term=terms[j],n_used=nrow(d),n_keep=keep,
   rand_HR_med=median(sims[,j],na.rm=T),rand_HR_q05=quantile(sims[,j],.05,na.rm=T),rand_HR_q95=quantile(sims[,j],.95,na.rm=T),rand_sig_share=mean(sims[,nt+j]<.05,na.rm=T))
}
res<-do.call(rbind,out); print(res,row.names=FALSE,digits=3); write.csv(res,'/private/tmp/incommon-randdrop.csv',row.names=FALSE)
