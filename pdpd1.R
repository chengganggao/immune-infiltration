see(unipd)
table(unipd$)
TNpd <- readcsv('pdTN.csv')
table(TNpd$pair)
scoreall <- readcsv('scoreall.csv')
scoreTN <- scoreall[rownames(TNpd),]

dat <- t(scoreTN)
rownames(dat) <- c("ABC Transporter","Amino Acid Metabolism","Aminosugar","Cholesterol", "Cofactor","Complex I","Complex IV","Detox","Fatty Acid","Glutathione","Glycan","Glycan Anchor","Glycan Degradation","Glycan Sulfate","Glycolysis","Heparin Sulfate","Hormone","Inositol Phosphate","Ion Transport","Krebs","Membrane Lipid","Multipurpose","NAD","Neurotransmitter","Nucleotide", "Other","Other Transport","Proton Transport","Purine","Pyrimidine","Redox","Signalling","Small Molecule Transport","Sphingolipid","Steroid","Sugar","Sulfate","Vitamin A")
group_list <- TNpd$pair
cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
n=t(scale(t(dat[cg,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n) 

pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')

pheatmap(n,cluster_cols = F,show_colnames =F,show_rownames = T,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')

library(KEGGREST)

exp <- readcsv('G:/Rstudio/project/pancreatic_cancer_infiltration/merge/expression/expression32.csv')
pdpd1 <- readcsv('pdpd2.csv')
exp2 <- t(exp)[rownames(pdpd1),]
see(exp2)
which(colnames(exp2)=='HLA-DRB4')
exp3 <- exp2[,c('HLA-DRB5','HLA-DRB4','HLA-DRB1','HLA-DRA','HLA-DQB1','HLA-DQA2','HLA-DQA1','HLA-DPB2','HLA-DPB1','HLA-DPA1','HLA-DOB','HLA-DOA','HLA-DMB','HLA-DMA')]
see(exp3)
colnames(exp3)
identical(rownames(exp3),rownames(pdpd1))
exp4 <- cbind(exp3,pdpd1)
colnames(exp4)[15:16]=c('os_status','os')
table(exp4$os_status)
exp4$os_status[which(exp4$os_status=='0')]='censored'
exp4$os_status[which(exp4$os_status=='1')]='uncensored'
write.csv(exp4,'expexp2.csv')
