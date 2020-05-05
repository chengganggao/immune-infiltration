rm(list = ls())
rm(gly)
library(vioplot)                                                    #引用包


pdn <- pdnt[which(pdnt$pair=='N'),]
pdt <- pdnt[which(pdnt$pair=='T'),]
rt <- rbind(pick(panall,pdn),pick(panall,pdt))
colnames(rt)




pdf("vioplot.pdf",height=8,width=50)              #保存图片的文件名称
par(las=1,mar=c(8,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))



path <- pdall[order(pdall$age),]
rt <- panall[rownames(path),]
rt <- rt[,1:38]

levels(pdall$age)
rt <- rbind(panall[which(pdall$mstage=='M0'),],panall[which(pdall$mstage=='M1'),])
table(pdall$mstage)
rt <- rt[,1:22]
see(rt)
normal=199                                                       #con样品数目
tumor=25
colnames(rt) <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
par(las=1,mar=c(8,6,3,3))
plot(x,y,
     xlim=c(0,64),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="infiltration",
     pch=21,
     col="white",
     xaxt="n")
#对每个免疫细胞循环，绘制vioplot
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
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.8,srt = 45,pos=2)
}



levels(pdall$mstage)
rt <- rbind(panall[which(pdall$mstage=='T1'),],panall[which(pdall$mstage=='T2'),],panall[which(pdall$mstage=='T1_or_T2'),],panall[which(pdall$mstage=='T3'),],panall[which(pdall$mstage=='T4'),])
table(pdall$mstage)
rt <- rt[,1:22]
see(rt)

T1=19                                                       #con样品数目
T2=86
T3=481
T4=21
colnames(rt) <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
par(las=1,mar=c(8,6,3,3))
plot(x,y,
     xlim=c(0,85),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="infiltration",
     pch=21,
     col="white",
     xaxt="n")

for(i in 1:ncol(rt)){
  if(sd(rt[1:G1,i])==0){
    rt[1,i]=0.001
  }
  if(sd(rt[(G1+1):(G1+G2),i])==0){
    rt[(G1+1),i]=0.001
  }
  if(sd(rt[(G1+G2):(G2+G3),i])==0){
    rt[(G1+G2),i]=0.001
  }
  G1Data=rt[1:G1,i]
  G2Data=rt[(G1+1):(G1+G2),i]
  G3Data=rt[(G1+G2):(G2+G3),i]
  vioplot(G1Data,at=4*(i-1),lty=1,add = T,col = 'seagreen3')
  vioplot(G2Data,at=4*(i-1)+1,lty=1,add = T,col = 'salmon1')
  vioplot(G3Data,at=4*(i-1)+2,lty=1,add = T,col = 'slateblue1')
  wilcoxTest=wilcox.test(G1Data,G2Data,paired = F)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(G1Data,G2Data))
  lines(c(x=4*(i-1)+0.2,x=4*(i-1)+0.8),c(mx,mx))
  text(x=4*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
  wilcoxTest=wilcox.test(G1Data,G3Data,paired = F)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(G2Data,G3Data))
  lines(c(x=4*(i-1)+1.2,x=4*(i-1)+1.8),c(mx-0.04,mx-0.04))
  text(x=4*(i-1)+1.5, y=mx-0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
  text(seq(1,85,4),-0.05,xpd = NA,labels=colnames(rt),cex = 0.8,srt = 45,pos=2)
}
dev.off()

