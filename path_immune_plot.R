library(vioplot) 

#2个变量
scorepdac <- pick(scorepdac,panpdac)
scorepdac <- scorepdac[rownames(panpdac),]
scorepanpdac <- cbind(panpdac[,1:22],scorepdac)
see(scorepanpdac)
colnames(rt)

rt <- scorepanpdac
rt <- rt[order(rt$Dendritic.cells.activated),]
rt <- rt[,23:60]
see(rt)

normal=474
tumor=474
colnames(rt) <- c("ABC Transporter","Amino Acid Metabolism","Aminosugar","Cholesterol", "Cofactor","Complex I","Complex IV","Detox","Fatty Acid","Glutathione","Glycan","Glycan Anchor","Glycan Degradation","Glycan Sulfate","Glycolysis","Heparin Sulfate","Hormone","Inositol Phosphate","Ion Transport","Krebs","Membrane Lipid","Multipurpose","NAD","Neurotransmitter","Nucleotide", "Other","Other Transport","Proton Transport","Purine","Pyrimidine","Redox","Signalling","Small Molecule Transport","Sphingolipid","Steroid","Sugar","Sulfate","Vitamin A")
par(las=1,mar=c(10,6,3,3))
x=NULL
y=NULL
plot(x,y,
     xlim=c(0,112),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Molecular metabolic pathway intensity",
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
  wilcoxTest=wilcox.test(normalData,tumorData,paired = T)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
  text(seq(1,112,3),-0.92,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
legend(27,0.96,legend = c('Dendritic cells activated (low) (215)','Dendritic cells activated (high) (215)'),col = c('seagreen3','salmon1'),pch = c(19,19),pt.cex = 2,bty = "n",cex = 0.8 )

