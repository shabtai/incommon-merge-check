library(survival)
f <- readRDS('/private/tmp/incommon-exact-kras-model.rds')
d <- readRDS('/private/tmp/incommon-exact-kras-input.rds')
src <- read.delim('/private/tmp/incommon-paper-sample.txt',comment.char='#',check.names=FALSE)
stopifnot(!anyDuplicated(src$SAMPLE_ID),!anyDuplicated(src$PATIENT_ID),!anyDuplicated(d$sample))
d$original_code <- src$ONCOTREE_CODE[match(d$sample,src$SAMPLE_ID)]
stopifnot(!anyNA(d$original_code),!anyDuplicated(src$PATIENT_ID[match(d$sample,src$SAMPLE_ID)]))
refit <- coxph(formula(f),data=d,ties=f$method)
stopifnot(isTRUE(all.equal(coef(refit),coef(f),tolerance=1e-10)),isTRUE(all.equal(vcov(refit),vcov(f),tolerance=1e-10)),isTRUE(all.equal(refit$loglik,f$loglik,tolerance=1e-10)),refit$n==f$n,refit$nevent==f$nevent)
cat('EXACT REPRODUCTION PASSED: coefficients, covariance, log-likelihood, N, deaths\n')
clean <- d[d$original_code=='PAAD',]
corrected <- coxph(formula(f),data=clean,ties=f$method)
cat('Original codes by class, pre missing-value exclusion:\n');print(table(d$original_code,d$class))
complete<-complete.cases(model.frame(formula(f),data=d,na.action=na.pass));cat('Original codes by class, analysed:\n');print(table(d$original_code[complete],d$class[complete]))
show <- function(f,label){s<-summary(f);data.frame(analysis=label,N=f$n,deaths=f$nevent,term=rownames(s$coef),HR=s$conf.int[,1],lower95=s$conf.int[,3],upper95=s$conf.int[,4],p=s$coef[,5])}
res<-rbind(show(refit,'published_merged'),show(corrected,'original_PAAD_only'))
print(res[grepl('class',res$term),],row.names=FALSE)
write.csv(res,'/private/tmp/incommon-exact-sensitivity.csv',row.names=FALSE)
write.csv(d[,c('sample','original_code','tumor_type','class','OS_MONTHS','OS_STATUS','age','TMB','type','FGA','sex')],'/private/tmp/incommon-exact-sensitivity-input.csv',row.names=FALSE)
