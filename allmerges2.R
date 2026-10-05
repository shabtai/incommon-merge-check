suppressMessages({library(survival)})
pub<-readRDS('/private/tmp/incommon-survival.rds'); pub<-pub[pub$fit_type=='general',]
pub<-pub[grepl('^(PAAD|LUAD|HCC|MEL|BLCA|OV|BRCA|UCEC|CRC)',pub$tumor_type),]
src<-read.delim('/private/tmp/incommon-paper-sample.txt',comment.char='#',check.names=FALSE)
named<-c(PAAD='PAAD',LUAD='LUAD',HCC='HCC',MEL='SKCM',BLCA='BLCA',OV='HGSOC',BRCA='IDC',UCEC='UEC',CRC='COAD')
g<-function(f,term){if(is.null(f)||!(term %in% names(coef(f)))||is.na(coef(f)[term]))return(c(NA,NA,NA,NA));s<-summary(f);c(s$conf.int[term,c(1,3,4)],s$coefficients[term,5])}
out<-list();comp<-list()
for(i in seq_len(nrow(pub))){p<-pub$cox_fit[[i]]; tt<-pub$tumor_type[i]; gene<-pub$gene[i]
 e<-environment(p$terms); d<-eval(p$call$data,envir=e)
 d$orig<-src$ONCOTREE_CODE[match(d$sample,src$SAMPLE_ID)]
 r<-coxph(formula(p),data=d,ties=p$method)
 ok<-isTRUE(all.equal(coef(r),coef(p),tolerance=1e-10)) && r$n==p$n && r$nevent==p$nevent
 used<-complete.cases(model.frame(formula(p),data=d,na.action=na.pass)); du<-d[used,]
 nm<-named[[sub(':.*','',tt)]]
 cm<-as.data.frame(table(orig=du$orig,class=du$class)); cm<-cm[cm$Freq>0,]; cm$gene<-gene; cm$tt<-tt; comp[[i]]<-cm
 r2<-tryCatch(coxph(formula(p),data=d[d$orig %in% nm,],ties=p$method),error=function(e)NULL)
 row<-data.frame(tt,gene,reproduced=ok,N=p$n,deaths=p$nevent,N_named=if(is.null(r2))NA else r2$n,named=nm,
   other_codes=paste(names(sort(table(du$orig[du$orig!=nm]),decreasing=TRUE)),collapse='+'),
   wt_n=sum(du$class=='WT'),wt_other=sum(du$class=='WT'&du$orig!=nm),mut_n=sum(du$class!='WT'),mut_other=sum(du$class!='WT'&du$orig!=nm))
 for(cl in c('PWK','BK','PMK')){a<-g(p,paste0('class',cl));b<-g(r2,paste0('class',cl))
  row[[paste0(cl,'_pub')]]<-a[1];row[[paste0(cl,'_pub_lo')]]<-a[2];row[[paste0(cl,'_pub_hi')]]<-a[3];row[[paste0(cl,'_pub_p')]]<-a[4]
  row[[paste0(cl,'_orig')]]<-b[1];row[[paste0(cl,'_orig_lo')]]<-b[2];row[[paste0(cl,'_orig_hi')]]<-b[3];row[[paste0(cl,'_orig_p')]]<-b[4]}
 out[[i]]<-row}
res<-do.call(rbind,out); write.csv(res,'/private/tmp/incommon-allmerges.csv',row.names=FALSE)
write.csv(do.call(rbind,comp),'/private/tmp/incommon-allmerges-composition.csv',row.names=FALSE)
cat('reproduced',sum(res$reproduced),'of',nrow(res),'\n')
