

#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.8")


setwd("C://Users//lexb4//Desktop//TCGAimmune//06.CIBERSORT")
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

panall <- readcsv('G:/Rstudio/project/pancreatic_cancer_infiltration/merge/expression/panall.csv')
panpdac <- pick(panall,pdpdac)
see(panpdac)
pansur <- pick(panpdac,pdsur)
pansur <- pansur[rownames(pdsur),]
colnames(pdsur)
pansur <- cbind(pdsur,pansur)
see(pansur)
see(scoresur)
see(pdpdac)
see(panpdac)
panall <- panall[which(panall$P.value<0.05),]
pdall <- pick(pdall,panall)
identical(rownames(pdpdac),rownames(panpdac))
panpdac <- panpdac[rownames(pdpdac),]
write.csv(pansur,'pansur.csv')
panpdac <- cbind(panpdac,pdpdac)
write.csv(pdall,'pdall1.csv')
pdpdac <- readcsv('TNpd.csv')




