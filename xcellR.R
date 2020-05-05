options(unzip = "internal")
devtools :: install_github(' dviraran / xCell ')


unipd <- readcsv('unipd.csv')
see(dat)
see(unipan)
panpdac <- dat
table(unipd$tissue)
ajpd <- unipd[which(unipd$tissue=='normal_PDAC'),]
normalpd <- unipd[which(unipd$tissue=='normal'),]
panaj <- pick(unipan,ajpd)
pannormal <- pick(unipan,normalpd)

colnames(pannormal)=c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")

par(mar=c(13,3,3,2))
boxplot(panpdac[,12:22],
        boxwex = 0.18, at = 1:11 - 0.22,pch=20,
        col = "salmon1",
        names = colnames(panpdac)[12:22],las=2,tck=0,
        main = "",
        xlab = "",
        ylab = "score",
        xlim = c(0.5, 11.5), ylim = c(-0.05, 0.72), yaxs = "i")
boxplot(panaj[,12:22], add = TRUE,
        boxwex = 0.18, at = 1:11 ,pch=20,
        names ='',
        col = 'mediumseagreen')
boxplot(pannormal[,12:22], add = TRUE,
        boxwex = 0.18, at = 1:11 + 0.22 ,pch=20,
        names ='',
        col = 'cadetblue1')


legend(9,0.7,legend = c('PDAC tissue','PDAC adjacent tissue','normal pancreas tissue'),col = c('salmon1','mediumseagreen','cadetblue1'),
       pch = c(19,19),pt.cex = 2,bty = "n",adj = )

dev.off()


rm(exp3)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

exp2 <- readcsv('G:/Rstudio/project/pancreatic_cancer_infiltration/merge/expression/expression3.csv')
exp2 <- t(exp2)
pan2 <- pick(panall,exp2)
pan2 <- pick(pan2,panpdac)
see(pan2)
pan_2 <- pan2
for (i in 1:22) {
  pan_2 <- pan_2[order(pan_2[,i]),]
  pan_2[,i][1:125]='low'
  pan_2[,i][126:250]='hTHEMIS'
}
dat <- t(as.matrix(pick(exp2,pan2)))
see(dat)
see(pan_2)
pan_2 <- pan_2[colnames(dat),]

dat=2^dat




group_list=group$group
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
exprSet=dat
contrast.matrix<-makeContrasts("hTHEMIS-low",
                               levels = design)

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

deg = deg(exprSet,design,contrast.matrix)
save(deg,file = 'deg.Rdata')
write.csv(deg,file = 'B.cells.memory_deg.csv')

#········································································
xzhdeg <- read.csv('GEO2R.csv')
xzhdeg <- xzhdeg[-which(duplicated(xzhdeg$Gene.symbol)),]
yellow <- read.csv('yellow.csv')
yellow <- yellow[-which(duplicated(yellow$symbol)),]
see(yellow)
yellow$data <- 1:177
namerow(xzhdeg)
deg <- xzhdeg
namerow(xzhdeg)
namerow(yellow)
rownames(yellow) <- yellow$symbol
see(yellow)
rownames(yellow)
rownames(xzhdeg)
rownames(xzhdeg) <- xzhdeg$Gene.symbol
deg <- pick(xzhdeg,yellow)
load(file = 'deg.Rdata')
head(deg)
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)




table(deg$g)
head(deg)
deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)


kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
upgene <- kk.up[1:3,1:2]
dotplot(kk.up );ggsave('kk.up.dotplot.png')
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.down)[,1:6]
downgene <- kk.down[1:3,1:2]
dotplot(kk.down );ggsave('kk.down.dotplot.png')
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]
dotplot(kk.diff );ggsave('kk.diff.dotplot.png')

kegg_diff_dt <- as.data.frame(kk.diff)
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
source('functions.R')
g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)

ggsave(g_kegg,filename = 'kegg_up_down.png',width = 12,heTHEMISt = 20,limitsize = F)

kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))

down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)

go_enrich_results <- lapply( g_list , function(gene) {
  lapply( c('BP','MF','CC') , function(ont) {
    cat(paste('Now process ',ont ))
    ego <- enrichGO(gene          = gene,
                    universe      = gene_all,
                    OrgDb         = org.Hs.eg.db,
                    ont           = ont ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.99,
                    qvalueCutoff  = 0.99,
                    readable      = TRUE)
    
    print( head(ego) )
    return(ego)
  })
})
save(go_enrich_results,file = 'go_enrich_results.Rdata')


load(file = 'go_enrich_results.Rdata')

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'/n'))
    png(fn,res=150,width = 1300)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}

install.packages('coplot')
BiocManager::install('coplot')

given.depth <- co.intervals(quakes$depth, number = 4, overlap = 0) #是将地震深度分为四个等级，类比肿瘤分期
coplot(lat ~ long | depth, data = quakes, given.v = given.depth, rows = 1, 
       panel = function(x,y,...) panel.smooth(x,y,span=0.7))

see(unipan)
colnames(unipan) <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils",'P value','Pearson Correlation','RMSE')
data1 <- unipan[,c(18,19,6,14,12,2)]
see(pan2)
see(dat)
dat1 <- cbind(pan2,t(dat))
identical(rownames(pan2),colnames(dat))
see(dat1)
colnames(dat1)[1:22]=c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
colnames(score_2)<- c("ABC Transporter","Amino Acid Metabolism","Aminosugar","Cholesterol", "Cofactor","Complex I","Complex IV","Detox","Fatty Acid","Glutathione","Glycan","Glycan Anchor","Glycan Degradation","Glycan Sulfate","Glycolysis","Heparin Sulfate","Hormone","Inositol Phosphate","Ion Transport","Krebs","Membrane Lipid","Multipurpose","NAD","Neurotransmitter","Nucleotide", "Other","Other Transport","Proton Transport","Purine","Pyrimidine","Redox","Signalling","Small Molecule Transport","Sphingolipid","Steroid","Sugar","Sulfate","Vitamin A")
dat1 <- cbind(dat1,score2)
score_2 <- score_2[rownames(pan2),]
see(score_2)

dat1$tissue <- 'PDAC'
immunegene <- readcsv('immunecell_gene.csv')
immgene <- immunegene[order(immunegene$T.cells.CD8),]
gene <- immgene[c(1:500,17467:17966),]
gene <- pick(gene,dat)
rownames(gene)


R=immunegene['THEMIS','T.cells.CD8']
P=immuneP['THEMIS','T.cells.CD8']
ymax=range(dat1$THEMIS)[2]
xmax=range(dat1$T.cells.CD8)[2]

ggplot(data=dat1,aes(x=T.cells.CD8,y=THEMIS))+geom_point()+geom_smooth()+geom_text(x = xmax-0.1, y = ymax, label = paste0('R=',R) ,parse = F)+geom_text(x = xmax-0.1, y = ymax-0.5, label = paste0('P=',P) ,parse = F)

library(fpc)
see(pansur)
cluster<-pansur[,26:48]
cluster <- cluster[which(cluster$P.value<0.05),]
ds<-dbscan(cluster,eps=0.42,MinPts = 5)
table(ds$cluster,iris$Species)

plot(ds,cluster)

library(fpc)
cluster.pamk<-pamk(cluster)

layout(matrix(c(1,2),1,2)) 
plot(cluster.pamk$pamobject)
layout(matrix(1))

library(cluster)
cluster.pamk<-clara(cluster,6)

clust <- cluster.pamk$clustering
clusterorder <- cluster.pamk$clustering
pansurcluster <- cbind(pansur0.05,clusterorder)
pansurcluster$clusterorder

cluster4 <- clusterorder


save(unipan,unipd,pdacpd,pansur,uniscore,file ='rowdata.Rdata')
load('rowdata.Rdata')


range(panpdac0.05$`T cells CD8`)
sum(panpdac0.05$`T cells CD8`<0.05)
