library(vioplot) 

pd <- readcsv('TNpd.csv')
see(pd)
xxx <- pick(unipd,rt)
table(xxx$tissue)
pan <- panall[rownames(pd),]

pan <- pan[rownames(pd),]
see(panpdac0.05)
identical(rownames(unipan),rownames(unipd))
scorepdac0.05 <- uniscore[rownames(panpdac0.05),]
see(scorepdac0.05)
panscorepdac0.05 <- cbind(panpdac0.05,scorepdac0.05)
rt <- panscorepdac0.05[order(panscorepdac0.05$Vitamin_A),][,1:22]

normal=474
tumor=474
colnames(rt) <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophpairs M0","Macrophpairs M1","Macrophpairs M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
par(las=1,mar=c(10,6,3,3))
x=NULL
y=NULL
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="infiltration",
     pch=21,
     col="white",
     xaxt="n")

for(i in 1:ncol(rt)){
  if(sd(rt[1:normal,i])==0){
    rt[1,i]=0.001
  }
  if(sd(rt[(normal+1):(normal+tumor),i])==0){
    rt[(normal+1),i]=0.001
  }
  normalData=rt[1:normal,i]
  tumorData=rt[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'seagreen3')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'salmon1')
  wilcoxTest=wilcox.test(normalData,tumorData,paired = F)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
legend(45,0.72,legend = c('Vitamin A level (low) (474)','Vitamin A level (high) (474)'),col = c('seagreen3','salmon1'),pch = c(19,19),pt.cex = 2,bty = "n",cex = 0.8 )

