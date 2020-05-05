rm(exp3)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

exp2 <- readcsv('G:/Rstudio/project/pancreatic_cancer_infiltration/merge/expression/expression3.csv')
exp2 <- t(exp2)
panpdac0.05 <- pick(unipan,pdacpd)
pan2 <- pick(panpdac0.05,exp2)
see(pan2)

pan_2 <- pan2
for (i in 1:22) {
  pan_2 <- pan_2[order(pan_2[,i]),]
  pan_2[,i][1:125]='low'
  pan_2[,i][126:250]='high'
}

score_2 <- score2
for (i in 1:38) {
  score_2 <- score_2[order(score_2[,i]),]
  score_2[,i][1:125]='low'
  score_2[,i][126:250]='high'
}
see(score_2)

dat <- t(as.matrix(pick(exp2,pan2)))
see(dat)
see(pan_2)
pan_2 <- pan_2[colnames(dat),]


table(pd2$prognosis)
group_list=pan_2$`T cells regulatory (Tregs)`
group_list=pd2$prognosis
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
exprSet=dat
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts("high-low",
                               levels = design)

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
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
write.csv(deg,file = 'T cells regulatory (Tregs)_deg.csv')

#········································································

load(file = 'deg.Rdata')
head(deg)
logFC_t=0
deg$g=ifelse(deg$adj.P.Val>0.05,'stable',
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
dotplot(kk.up );ggsave('kk.up.dotplot.png')
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.down)[,1:6]
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

ggsave(g_kegg,filename = 'kegg_up_down.png',width = 12,height = 10,limitsize = F)

kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
kk_gsea=gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
ggsave(kk_gsea,filename = 'kk_gsea.png')

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
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1300)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}

