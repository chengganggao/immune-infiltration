library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) 
library(GSEABase)
rm()
d='G:/Rstudio/project/pancreatic_cancer_infiltration/gmt'
gmts=list.files(d,pattern = 'all')
gmts <- gmts[1]
geneset <- getGmt(file.path(d,gmts))  

exppathwayscore <- lapply(expall,function(x)gsva(as.matrix(x), geneset, 
                 mx.diff=FALSE, verbose=FALSE, 
                 parallel.sz=1))
exppathwayscore <- exppatnwayscore
see(exppathwayscore[[1]])
rownames(exppathwayscore[[1]])
exppathway <- exppathwayscore[[1]]
for (i in 2:30) {
  i
  exppathway <- cbind(exppathway,exppathwayscore[[i]])
}
see(exppathway)

dim(scoresur)
scoresur <- pick(t(exppathway),pansur)
scoresur <- scoresur[rownames(pansur),]
scoresur <- cbind(scoresur,pansur[,26:50])
write.csv(scoresur,'scoresur.csv')
colnames(scoresur)


-





























if(T){
    es_max <- lapply(gmts, function(gmtfile){ 
      #gmtfile=gmts[8];gmtfile
      geneset <- getGmt(file.path(d,gmtfile))  
      es.max <- gsva(X, geneset, 
                     mx.diff=FALSE, verbose=FALSE, 
                     parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max, function(es.max){
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                                     levels = design)
      contrast.matrix<-makeContrasts("Tumor-Normal",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
  } 
  
  gmts
  
  save(es_max,es_deg,file='gsva_msigdb.Rdata')
  
  load(file='gsva_msigdb.Rdata')
  
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=2
    print(i)
    dat=es_max[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.01 & abs(df$logFC) > 0.3,]
    print(dim(df))
    if(nrow(df)>5){
      n=rownames(df)
      dat=dat[match(n,rownames(dat)),]
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(dat)
      rownames(dat)=substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
      
    }
  })
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max)
  df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df,file = 'GSVA_DEG.csv')


