library(ggplot2)
ggplot(mpg,aes(x=class)) + geom_bar() 
ggplot(mpg,aes(x=class,y=displ)) + geom_bar(stat="identity") 
ggplot(mpg,aes(x=class)) + geom_bar(aes(weight=displ)) 
ggplot(mpg,aes(x=class)) + geom_bar(aes(weight=rep(1,length(class)))) 

ggplot(mpg,aes(x=class)) + geom_bar(aes(fill=factor(cyl)),position="stack") 
ggplot(mpg,aes(x=class)) + geom_bar(aes(fill=factor(cyl)),position="dodge") 

ggplot(mpg,aes(x=class)) + geom_bar(aes(fill=factor(cyl)),position="fill") 

see(panall2)
panall <- readcsv('G:/Rstudio/project/pancreatic_cancer_infiltration/merge/expression/panall.csv')
unipan <- pick(panall,unipd)
unipan <- unipan[rownames(unipd),]
see(unipan)
see(unipd)
unipanbar <- cbind(unipan,unipd$batch)
table(unipd$batch)
write.csv(unipanbar,'unipanbar.csv')
see(unipanbar)
unipanbar$p <- 3
unipanbar$p[which(unipanbar$P)]=1
colnames(unipanbar)
unipb <- unipanbar[,26:27]
number <- c(1:4,6:22,24:30)
vector=sapply(number,function(x)sum(unipb$p==1 & unipb$batch==x)/ sum(unipb$batch==x))
studyname <- studyname[order(vector)]
tnum <- seq(0.01,0.28,by=0.01)
for (i in 1:28) {
  unipb[which(unipb$batch==num[i]),1] <- tnum[i]
}
unipb <- unipb[order(unipb$batch),]
for (i in 1:28) {
  unipb[which(unipb$batch==tnum[i]),1] <- studyname[i]
}
table(unipb$batch)
studyname2 <- studyname[c(1,2,4,5,3,6:28)]
unipb$batch <- factor(unipb$batch, levels = studyname)
unipb[which(unipb$p==3),2]='p≥0.05'

g=ggplot(unipb,aes(x=batch)) + geom_bar(aes(fill=factor(p)),position="fill") 
g=g+theme(axis.text.x=element_text(angle=90,size=8))
g=g+ylab('pattern of p-value')
g

unipanbar$batch
colnames(unipanbar)[26] <- 'batch'
table(unipanbar$batch)
df <- data.frame('immunecell'=colnames(unipanbar)
                   )
for (i in c(1:4,6:22,24:30)) {
  df1 <- unipanbar[which(unipanbar$batch==i),]
  mean1 <- apply(t(df1), 1, mean)
  df <- cbind(df,mean1)
}
colnames(df)[2:29]=c(1:4,6:22,24:30)
df <- namerow(df)
df <- t(df)

df <- as.data.frame(df)
class(df)
df <- df[1:22,]
df$immunecell=c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
library(reshape2)
head(mpg)
View(mpg)
df <- df[,c(1:22,27)]
df_long <- reshape(df,varying = list(2:29),v.names = 'immunecell',idvar = 'cyl',direction = 'long')
see(df)
see(df_long)
colnames(df_long)[2]='class'
study <- read.csv('study.csv')
study$name<- paste0(study$Accession.number.Source,'(N=',study$Number.of.patients,')')

write.csv(df_long,'df_long.csv')
immunename <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
studyname <- study$name

studyname <- studyname[order(vector)]

options(stringsAsFactors = FALSE) 
study1 <- as.data.frame(t(readcsv('study1.csv')))
study1$class <- study$name
colnames(study1)[1:2] <- c('count','survival')
study1$count <- as.numeric(study1$count)
study1$class <- factor(study1$class, levels = studyname)
study1$survival[which(study1$survival=='no')]=0
study1$survival <- as.numeric(study1$survival)
g=ggplot(study1,aes(x=class,y = survival,)) + 
  geom_bar(stat = "identity") 
g=g+theme(axis.text.x=element_text(angle=90,size=8))
g+ylab('survival data available')

head(mtcars)

plot(1:29,xlab=studyname)
see(pansurcluster)
clusterbar <- pansurcluster[,c(26:47,51)]
clusterbar$clusterorder <- cluster5
df <- NA
table(clusterbar$clusterorder)
for (i in c(1:5)) {
  df1 <- clusterbar[which(clusterbar$clusterorder==i),]
  mean1 <- apply(t(df1), 1, mean)
  df <- cbind(df,mean1)
}
df <- as.data.frame(df[-23,2:6])
colnames(df)=c("cluster1","cluster2",'cluster3','cluster4','cluster5')
df$cell <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")

clusterbar1 <- as.data.frame(t(clusterbar[which(clusterbar$clusterorder==5),-23]))
clusterbar1$cell <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
see(clusterbar1)
clusterbar1_long <- reshape(clusterbar1,varying = 1:ncol(clusterbar1)-1,v.names = 'immunecell',idvar = 'cell',direction = 'long')
g=ggplot(clusterbar1_long,aes(x=time)) + geom_bar(aes(weight=immunecell,fill=factor(cell)),width = 1,position="fill")
#g=g+scale_fill_gradientn(colors=c('antiquewhite4','aquamarine','bisque2','blueviolet','brown','burlywood1','chartreuse1','chocolate1','cornflowerblue','cyan4','slategray4','darkgray','deeppink','gold','greenyellow','indianred1','lightskyblue2','orangered','royalblue4','tan1','darkturquoise','cadetblue'))
g+ylab('pattern of immune cell proportion')

df_long <- reshape(df,varying = 1:5,v.names = 'immunecell',idvar = 'cell',direction = 'long')

ggplot(df_long,aes(x=time)) + geom_bar(aes(weight=immunecell,fill=factor(cell)),position="fill") 
g=g+theme(axis.text.x=element_text(angle=90,size=8))
g+ylab('pattern of immune cell proportion')



