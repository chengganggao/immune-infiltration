library(KEGGREST)

upgene
downgene
gs<-keggGet('hsa04659 ')
names(gs[[1]])
gs[[1]]$GENE
class(gs[[1]]$GENE)
strsplit(gs[[1]]$GENE,';')
class(strsplit(gs[[1]]$GENE,';'))
genes <- unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';')[[1]][1]))
genes <-  genes[1:length(genes)%%2 ==0]

geneset6 <- data.frame(
  'Th17 cell differentiation'=genes
)

library("plyr")
list1 <- list(geneset,geneset2,geneset3,geneset4,geneset5,geneset6)
Macrophages.M0_kegggene=do.call(rbind.fill,list1)
Macrophages.M0_kegggene <- merge(geneset,geneset2,)









