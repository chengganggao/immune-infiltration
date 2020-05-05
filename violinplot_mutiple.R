library(vioplot) 

rt <- rbind(panall[which(pdall$grade=='G1'& pdall$tissue=='normal_PDAC'),],panall[which(pdall$grade=='G2'& pdall$tissue=='normal_PDAC'),],panall[which(pdall$grade=='G3'& pdall$tissue=='normal_PDAC'),])

see(scorepdac0.05)
see(clusterbar)

scorebar <- cbind(scorepdac0.05[rownames(clusterbar),],clusterbar[,23])
table(clusterbar$clusterorder)
rt <- rbind(clusterbar[which(clusterbar$clusterorder=='1'),],clusterbar[which(clusterbar$clusterorder=='2'),],clusterbar[which(clusterbar$clusterorder=='3'),],clusterbar[which(clusterbar$clusterorder=='4'),],clusterbar[which(clusterbar$clusterorder=='5'),])
rt <- rt[,1:22]
see(rt)
lev1=54                                                      
lev2=144
lev3=82
lev4=108
lev5=94
colnames(rt) <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
par(las=1,mar=c(8,6,3,3))
x=NULL
y=NULL
plot(x,y,
xlim=c(0,130),ylim=c(min(rt)-0.05,max(rt)+0.02),
main="",xlab="", ylab="infiltration",
pch=21,
col="white",
xaxt="n")
for(i in 1:ncol(rt)){
 if(sd(rt[1:lev1,i])==0){
    rt[1,i]=0.001
 }
 if(sd(rt[(lev1+1):(lev1+lev2),i])==0){
    rt[(lev1+1),i]=0.001
 }
 if(sd(rt[(lev1+lev2):(lev2+lev3),i])==0){
    rt[(lev1+lev2),i]=0.001
 }
 if(sd(rt[(lev2+lev3):(lev3+lev4),i])==0){
        rt[(lev2+lev3),i]=0.001
 }
 if(sd(rt[(lev3+lev4):(lev4+lev5),i])==0){
        rt[(lev3+lev4),i]=0.001
 }
 
 lev1Data=rt[1:lev1,i]
 lev2Data=rt[(lev1+1):(lev1+lev2),i]
 lev3Data=rt[(lev1+lev2):(lev2+lev3),i]
 lev4Data=rt[(lev2+lev3):(lev3+lev4),i]
 lev5Data=rt[(lev3+lev4):(lev4+lev5),i]
 
 vioplot(lev1Data,at=6*(i-1),lty=1,add = T,col = 'seagreen3')
 vioplot(lev2Data,at=6*(i-1)+1,lty=1,add = T,col = 'salmon1')
 vioplot(lev3Data,at=6*(i-1)+2,lty=1,add = T,col = 'slateblue1')
 vioplot(lev4Data,at=6*(i-1)+3,lty=1,add = T,col = 'antiquewhite4')
 vioplot(lev5Data,at=6*(i-1)+4,lty=1,add = T,col = 'gold')
 
 wilcoxTest=wilcox.test(lev1Data,lev2Data,paired = F)
 p=round(wilcoxTest$p.value,3)
 mx=max(c(lev1Data,lev2Data))
 lines(c(x=6*(i-1)+0.2,x=6*(i-1)+0.8),c(mx,mx))
 text(x=6*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
 wilcoxTest=wilcox.test(lev1Data,lev3Data,paired = F)
 p=round(wilcoxTest$p.value,3)
 mx=max(c(lev2Data,lev3Data))
 lines(c(x=6*(i-1)+1.2,x=6*(i-1)+1.8),c(mx-0.04,mx-0.04))
 text(x=6*(i-1)+1.5, y=mx-0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
 wilcoxTest=wilcox.test(lev1Data,lev4Data,paired = F)
 p=round(wilcoxTest$p.value,3)
 mx=max(c(lev2Data,lev3Data))
 lines(c(x=6*(i-1)+2.2,x=6*(i-1)+2.8),c(mx-0.08,mx-0.08))
 text(x=6*(i-1)+2.5, y=mx-0.06, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
 wilcoxTest=wilcox.test(lev1Data,lev5Data,paired = F)
 p=round(wilcoxTest$p.value,3)
 mx=max(c(lev2Data,lev3Data))
 lines(c(x=6*(i-1)+3.2,x=6*(i-1)+3.8),c(mx-0.1,mx-0.1))
 text(x=6*(i-1)+3.5, y=mx-0.08, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
 
 text(seq(3,130,6),-0.1,xpd = NA,labels=colnames(rt),cex = 0.8,srt = 45,pos=2)
}
legend(100,0.75,legend = c('Cluster5 (94)','Cluster1 (54)','Cluster2 (144)','Cluster3 (82)','Clster4 (108)'),col = c('seagreen3','salmon1','slateblue1','antiquewhite4','gold'),
       pch = c(19,19),pt.cex = 1.5,bty = "n",cex = 0.8 )
