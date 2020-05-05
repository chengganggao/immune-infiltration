#正态性及方差分析
rm(list = ls()) 
load(file = 'step1-output.Rdata')

library(car)

data <- t(dat)[,-22]
data <- data[,-28]
data <- data[,-39]
colnames(data)
colnames(tissueexp)
tissueexp <- data.frame(phe$tissue,stringsAsFactors = T)
colnames(tissueexp)[1] <- 'tissue'
#tissuedif <- data.frame('diff'=1,'lwr'=1,'upr'=1,'p_adj'=1,'cell'=1)
tissueexp$tissue <- as.factor(tissueexp$tissue)
nlevels(tissueexp$tissue)
i <- 1
for (i in 2:64) {
  n <- as.character(powerTransform(data[,i]+0.0001))
  v=as.numeric(substr(n[7],regexpr("=",n[7])+1,regexpr("=",n[7])+11))
  #tempname <- (data[,i]+0.0001)^v
  lam <- c(lam,v)
}
write.csv(tissueexp,file = 'tissueexp.csv')

for (i in 1:66) {
  n <- as.character(powerTransform(data1[,i]+0.0001))
  v=as.numeric(substr(n[7],regexpr("=",n[7])+1,regexpr("=",n[7])+11))
  tempname <- (data1[,i]+0.0001)^v
  tissueexp <- cbind(tissueexp,tempname)
  fit <- aov(tempname~tissue,data=tissueexp, var.equal = T)
  c <- as.data.frame(TukeyHSD(fit)$tissue)
  colnames(c)[4] <- 'p_adj'
  if (sum(c[,4]<0.05)>0) {
    tempdat <- cbind(c[which(c[,4]<0.05),],'cell'=rep(colnames(data1)[i], times=as.numeric(sum(c[,4]<0.05))))
    tissuedif <- rbind(tissuedif,tempdat)
  }else NULL
  colnames(tissueexp)[i+1] <- colnames(data1)[i]
}
tissuedif <- tissuedif[-1,]
write.csv(tissuedif,file = 'tissuedif.csv')



#克利夫兰条件图：
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
dat <- t(as.data.frame(dat))
dat[1:4,1:4]
head(m)
m <- data.frame(phe$days_to_death,phe$stage)
data_stage <- cbind(dat,m)
head(data_stage)
colnames(data_stage)
colnames(data_stage)[68:69] <- c('days','stage')

library(coplot)

data9 <- as.data.frame(cbind(t(es_max),data8[,c(65,68)]))
colnames(data9)[39] <- 'ImmuneScore'
data8$tissue <- as.factor(data8$tissue)
data9 <- as.data.frame(cbind(data8,datag))
head(data6)
i <- 1
for (i in 1:38) {
  name <- paste(colnames(data9)[i],'.png')
  png(filename = name)
  coplot(ImmuneScore ~ data9[,i] | tissue, data = data9, rows = 1,
  panel = function(x,y,...)panel.smooth(x,y,span=1,col = rainbow(20),
                                             pch = 19,col.smooth = 'brown'))
 dev.off()
  }

class(data7)

dev.off()

#上下柱状图
Ncell3 <- Ncell1[order(Ncell1)]
Ncell4 <- Ncell2[order(Ncell1)]
Nname1 <- Nname[order(Ncell1)]
Nsd3<- Nsd1[order(Ncell1)]
Nsd4 <- Nsd2[order(Ncell1)]



par(mar=c(9,3,3,2))
h1=barplot(Ncell3,col = 'salmon1',ylim = c(0,0.7),axes = T,)
axis(side = 1,at = h1,labels = Nname1,tick = F,las=3)
h2 = barplot(Ncell4,col = 'mediumseagreen',add = T)

legend(0.5,0.6,legend = c('Normal','Tumor'),col = c('salmon1','mediumseagreen'),
       pch = c(19,19),pt.cex = 2,bty = "n",adj = )
segments(h1,Ncell3-Nsd3/5,h1,Ncell3+Nsd3/5,lty = 1.2,col = 'orangered')
arrows(h1,Ncell3-Nsd3/5,h1,Ncell3+Nsd3/5,code = 3,length = 0.1,angle = 90,col = 'orangered')
segments(h1,Ncell4-Nsd4/5,h1,Ncell4+Nsd4/5,lty = 1.2,col = 'darkgreen')
arrows(h1,Ncell4-Nsd4/5,h1,Ncell4+Nsd4/5,code = 3,length = 0.1,angle = 90,col = 'darkgreen')


Tcell3 <- Tcell1[order(Tcell1)]
Tcell4 <- Tcell2[order(Tcell1)]
Tname1 <- Tname[order(Tcell1)]
Tname2 <- Tname1[seq(17,1,-1)]
Tsd3<- Tsd1[order(Tcell1)]
Tsd4 <- Tsd2[order(Tcell1)]



par(mar=c(9.5,3,3,2))
h1=barplot(Tcell3,col = 'mediumseagreen',ylim = c(0,0.7),axes = T,)
axis(side = 1,at = h1,labels = Tname2,tick = F,las=2)
h2 = barplot(Tcell4,col = 'salmon1',add = T)

legend(0.5,0.6,legend = c('Normal','Tumor'),col = c('salmon1','mediumseagreen'),
       pch = c(19,19),pt.cex = 2,bty = "n",adj = )
segments(h1,Tcell3-Tsd3/5,h1,Tcell3+Tsd3/5,lty = 1.2,col = 'darkgreen')
arrows(h1,Tcell3-Tsd3/5,h1,Tcell3+Tsd3/5,code = 3,length = 0.1,angle = 90,col ='darkgreen')
segments(h1,Tcell4-Tsd4/5,T1,Tcell4+Tsd4/5,lty = 1.2,col =  'orangered')
arrows(h1,Tcell4-Tsd4/5,h1,Tcell4+Tsd4/5,code = 3,length = 0.1,angle = 90,col =  'orangered')


#t.test
i <- 1
for (i in 2:64) {
  name_i=colnames(tissueexp)[i+1]
  tissueexp2=tissueexp[,c(1,i+1)]
  colnames(tissueexp2)=c('tissue','gene')
  a=t.test(gene~tissue,tissueexp2,paired=T)
  print(name_i)
  m <- c(m,a$p.value)
}
View(tissueexp)
length(tissueexp$tissue)
as.character(a)

data3 <- data[,which(m<0.05)]
lam1 <- lam[which(m<0.05)]
dim(data2)
b <- as.data.frame(tissueexp[,1])
class(b)
data4 <- cbind(data3,b)
dim(data3)
data3 <- (data3)^0.25
colnames(data4)[37] <- 'tissue'
View(data4)
data3T <- data3[seq(from = 1, to = 89, by = 2),]
data3N <- data3[seq(from = 2, to = 90, by = 2),]
mean1 <- aggregate(.~tissue,data = data4,mean)
se1 <- aggregate(.~tissue, data = data4,se)

for (i in 2:35) {
  if(mean1[1,i]>mean1[2,i]) {
    b1 <- mean1[1,i]
    b2 <- mean1[2,i]
    b3 <- colnames(mean1)[i]
    b4 <- sd1[1,i]
    b5 <- sd1[2,i]
    Ncell1 <- c(Ncell1,b1)
    Ncell2 <- c(Ncell2,b2)
    Nname <- c(Nname,b3)
    Nsd1 <- c(Nsd1,b4)
    Nsd2 <- c(Nsd2,b5)
  } else {
  c1 <- mean1[2,i]
  c2 <- mean1[1,i]
  c3 <- colnames(mean1)[i]
  c4 <- sd1[2,i]
  c5 <- sd1[1,i]
  Tcell1 <- c(Tcell1,c1)
  Tcell2 <- c(Tcell2,c2)
  Tname <- c(Tname,c3)
  Tsd1 <- c(Tsd1,c4)
  Tsd2 <- c(Tsd2,c5)
  }
}
Ncell1 <- NULL
Ncell2 <- NULL
Tcell1 <- NULL
Tcell2 <- NULL
Nname <- NULL
Tname <- NULL
Nsd1 <- NULL 
Nsd2<- NULL
Tsd1 <- NULL
Tsd2<- NULL

barplot(-log(Ncell),space = 0,col = rainbow(18))
dev.off()
class(mean1[1,-1])
n <- as.matrix(mean1[1,-1])
head(n)               
View(tissueexp)
View(data1)


#颜色选取

install.packages('RClorBrewer')
BiocManager::install('RClorBrewer')

library(RColorBrewer)
display.brewer.all(type='seq')
age <- rnorm(100,45,10)
hist(age,col = brewer.pal(9,'BuGn'))
display.brewer.all(type = 'div')
age <- rnorm(100,45,10)
hist(age,col = brewer.pal(11,'BrBG'))
display.brewer.all(type = 'qual')
value <- sample(20:100,10)
barplot(value,col = brewer.pal(10,'Set3'))
devtools::install_github("daattali/colourpicker")
CPCOLS <- "#006400"
plot(1:5,col = CPCOLS)


# 箱盒图
es_deg <- read.csv('es_deg.csv',header = T)
rownames(es_deg) <- es_deg$X
es_deg <- es_deg[,-1]
es_max <- read.csv('gsvascore.csv',header = T)
rownames(es_max) <- es_max$X
es_max <- es_max[,-1]

datag <- as.data.frame(t(es_max))
data5 <- datag[,rownames(es_deg)[c(1:7,12,15,17,18)]]
data6 <- as.data.frame(cbind(data5,phe$tissue))
colnames(data6)[12] <- 'tissue'

library(reshape2)
data6_m <- melt(data = exp, id.vars = "tissue")
data6_m <- reshape(data6, idvar = "tissue",varying = ,
          v.names = "score", direction = "long")
datag <- as.data.frame(es.max)

class(datag)
range(datag)
k1 <- 

par(mar=c(13,3,3,2))
boxplot(value ~ variable, data = data6_m,
        boxwex = 0.25, at = 1:11 - 0.2,pch=20,
        #at 参数定义了盒形图中盒子的位置，传入一个向量，所以，
        #最终三个盒子的横坐标分别是 0.8，1.8，2.8
        subset = tissue == "N", col = "salmon1",
        names = colnames(data6)[1:11],las=3,tck=0,
        main = "metabolic pathway",
        xlab = "",
        ylab = "score",
        xlim = c(0.5, 11.5), ylim = c(-0.8, 0.8), yaxs = "i")
boxplot(value ~ variable, data = data6_m, add = TRUE,
        #在当前的图形上，添加一个新的盒形图
        boxwex = 0.25, at = 1:11 + 0.2,pch=20,
        names ='',
        subset = tissue == "T", col = 'mediumseagreen')
legend(0.4,-0.5,legend = c('Normal','Tumor'),col = c('salmon1','mediumseagreen'),
       pch = c(19,19),pt.cex = 2,bty = "n",adj = )

dev.off()
