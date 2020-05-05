rm(UniVar2)

library(survival)
library(plyr)

coxdata <- readcsv('cox.csv')
pansur <- pansur[rownames(coxdata),]
pan_sur <- pansur
pan_sur[which(pan_sur$P.value>=0.05),]=NA
colnames(pan_sur)
coxdata <- cbind(pan_sur[,26:47],coxdata)
coxdata <- pick(pansur0.05, dat)

see(coxdata)
colnames(coxdata)
coxdata <- coxdata[,c(13,14,26:47)]
coxdata <- read.csv('Multi.csv')
coxdata <- readcsv('pansur0.05.csv')
coxdata <- readcsv('multidata1.csv')
sur <- Surv(time = coxdata$os,event = coxdata$os_status)

UniCox <- function(X){
  FML <- as.formula(paste0('sur~',X))
  tempcox <- coxph(FML,data = coxdata)
  tempsum <- summary(tempcox)
  HR <- round(tempsum$coefficients[,2],2)
  pValue <- round(tempsum$coefficients[,5],3)
  CI <- paste0(round(tempsum$conf.int[,3:4],2),collapse = '-')
  Unicox <- data.frame('Characteristics'=X,
                       'Hazard Ratio'=HR,
                       'CI 95'=CI,
                       'p Value'=pValue)
  return(Unicox)
}

UniCox <- function(X){
  FML <- as.formula(paste0('sur~',X))
  tempcox <- coxph(FML,data = coxdata)
  tempsum <- summary(tempcox)
  HR <- round(tempsum$coefficients[,2],2)
  pValue <- round(tempsum$coefficients[,5],3)
  CI <- paste0(round(tempsum$conf.int[,3:4],2),collapse = '-')
  N <- tempsum$n
  Events <- tempsum$nevent
  Unicox <- data.frame('Characteristics'=X,
                       'Hazard Ratio'=HR,
                       'CI 95'=CI,
                       'p Value'=pValue,
                       "N"=N,
                       'Events'=Events)
  return(Unicox)
}


VarNames <- colnames(coxdata)[1:13]
UniVar <- lapply(VarNames,UniCox)
UniVar <- ldply(UniVar,data.frame)
UniVar
write.csv(UniVar,'UniVar_1.csv')

MultiCoxdata <- coxdata[,UniVar$Characteristics[UniVar$p.Value<0.05]]
see(MultiCoxdata)
library(mice)
md.pattern(MultiCoxdata)
aggr_plot <- aggr(MultiCoxdata, col = c('red', 'green'), numbers = TRUE, sortVars = TRUE,labels= colnames(MultiCoxdata), cex.axis = 0.3,gap = 3)

MultiCoxdata <- cbind(MultiCoxdata,coxdata$os)
colnames(MultiCoxdata)[14]='os'

multicox
mice_mod <- mice(MultiCoxdata[7:14], method = 'rf',family= 'logistic')
mice_output <- complete(mice_mod)
pdata <- read.csv('MultiCoxdata.csv')

MultiCoxdata <- cbind(MultiCoxdata[,1:6],pdata[,2:10])
scoresur <- pick(scoreall,MultiCoxdata)
scoresur <- scoresur[rownames(MultiCoxdata),]
colnames(scoresur)
MultiCoxdata1 <- cbind(MultiCoxdata,scoresur[,c(22,35,8,18,2,5,30,15)])
see(MultiCoxdata1)
colnames(MultiCoxdata1)
MultiCoxdata1 <- readcsv('Multi.csv')
UniVar=readcsv('UniVar.csv')
table(MultiCoxdata1$os_status)
sur <- Surv(time = MultiCoxdata1$os,event = MultiCoxdata1$os_status)
fml <- as.formula(paste0('sur~',paste0(colnames(MultiCoxdata1)[1:12],collapse='+')))
MultiCox <- coxph(fml,data = MultiCoxdata1)
Multisum<- summary(MultiCox)

MultiName <- as.character(colnames(MultiCoxdata1)[1:12])
MHR <- round(Multisum$coefficients[,2],2)
MpValue <- round(Multisum$coefficients[,5],3)
MCIL <- round(Multisum$conf.int[,3],2)
MCIU <- round(Multisum$conf.int[,4],2)
MCI <- paste0(MCIL,'-',MCIU)
Multicox <- data.frame('Characteristics'=MultiName,
                     'Hazard Ratio'=MHR,
                     'CI 95'=MCI,
                     'p Value'=MpValue)
Multicox
write.csv(Multicox,'Multicoxfinal.csv')
Final <- merge.data.frame(UniVar,Multicox,by='Characteristics',
                          all=T,sort=T)
write.csv(Final,'finalcox.csv')

library(party)
write.csv(MultiCoxdata1,'Multi.csv')
control <- ctree_control(teststat = c("quad", "max"), 
                          testtype = c("Bonferroni", "MonteCarlo", 
                                       "Univariate", "Teststatistic"), 
                          mincriterion = 0.7, minsplit = 20, minbucket = 7, 
                          stump = FALSE, nresample = 9999, maxsurrogate = 0, 
                          mtry = 0, savesplitstats = TRUE, maxdepth = 0, remove_weights = FALSE)
tree <- ctree(Surv(os,os_status)~B.cells.memory+NK.cells.activated+Dendritic.cells.activated+grade+nstage+resection_margin+cdkn2amut ,data = MultiCoxdata1,controls = control)
plot(tree)






