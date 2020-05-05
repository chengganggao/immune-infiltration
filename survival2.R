library(forestplot)
library(tidyverse)
library(survival)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

source("function_subgroup_survival_analysis.R")
{cell <- c("B.cells.naive","B.cells.memory","Plasma.cells","T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting","NK.cells.activated","Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated","Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils")

input <- pansur
see(input)
input$pro[1:794]='ok'
input <- input[which(input$P.value<0.05),]


data <- lapply(cell,function(x)Subgroup_survival_analysis(pdata = input,
                                    time ="os", status = "os_status",
                                    variable = "pro",object = x))
foredata <- data[[1]][1,]
for (i in 2:22) {
  foredata <- rbind(foredata,data[[i]][1,])
  }
rownames(foredata) <- colnames(input)[26:47]
foredata

rownames(foredata) <- c("B cells naive","B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta","NK cells resting","NK cells activated","Monocytes","Macrophages M0","Macrophages M1","Macrophages M2","Dendritic cells resting","Dendritic cells activated","Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")
foredat <- as.matrix(foredata[which(foredata$CI_up_0.95<10000),])
foredat <- cbind(as.data.frame(foredat),mean=apply(foredat[,3:4], 1,mean))
foredat <- foredat[order(foredat$mean),]
foredat_inf <- foredata[c(which(foredata$CI_up_0.95>10000),which(foredata$CI_low_0.95>10000)),]
foredat_inf <- foredat_inf[1:9,]
foredat_inf <- foredat_inf[order(foredat_inf$CI_low_0.95),]
foredat_inf <- cbind(foredat_inf,mean=rep(1000,nrow(foredat_inf)))
text <- rbind(foredat,foredat_inf)
colnames(text)[3:4]=c('lower','upper')
cochrane_from_rmeta <- rbind(c(NA,NA,NA),text[,c(5,3,4)])
cochrane_from_rmeta[(nrow(foredat)+2):23,3] <- rep(2000,nrow(foredat_inf))
cochrane_from_rmeta[nrow(cochrane_from_rmeta),2]=0

tabletext<-cbind(c("Cell type",rownames(text)),
                 c('OR',text$HR),
                 c('95% CI',paste0('(',text$lower,',',text$upper,')'))
)

{as <- readcsv('as.csv')
as<- cbind(as,mean=apply(as, 1,mean))
as <- as[,c(3,1,2)]
as <- rbind('1'=NA,
            as)
as <- as[order(as$mean),]
cochrane_from_rmeta <- as

bs <- UniVar[23:35,c(1:3,5:6)]
rownames(bs) <- bs[,1]
bs <- bs[rownames(as),]
class(bs)
bs[,1] <- c("1"  ,                 "Perineural invasion", "Gender"       ,       "Age"    ,            
"KRAS mutation"       ,      "Location"        ,    "Tumor size"       ,       "TP53 mutation"   ,         
"Node status"    ,          "Grade"       ,        "Matastasis"     ,         "Resection margin"   ,
"CDKN2A mutation"       ,    "Vessel invasion" )
bs[1,] <- c('Characteristics','HR',"95%CI",'N','Events')
bs <- cbind(bs,
            'Categories'=c('Categories','Negative, Positive','male,female','<=65y,>65y','Negative, Positive','head,body/tail','T1,T2,T3,T4','Negative, Positive','Negative, Positive','G1,G2,G3','Negative, Positive','Negative, Positive','Negative, Positive','Negative, Positive'))
tabletext <- bs}
as <- readcsv('as1.csv')
bs <- readcsv('bs1.csv')
bs <- bs[rownames(as),]
bs <- cbind('Variables'=rownames(bs),bs)
bs[1,] <- c('Variables','HR',"95%CI")
cochrane_from_rmeta <- as
tabletext <- bs

coxdata <- readcsv('Multicoxfinal.csv')
colnames(coxdata)
cochrane_from_rmeta <- coxdata[,5:7]
tabletext <- coxdata[,1:4]

forestplot(tabletext, 
           graph.pos = 2,
           hrzl_lines = gpar(col="#444444"),
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,50),TRUE),
           clip=c(0.05,20), 
           xlog=TRUE,
           title = 'Prognostic associations of multifactors',
           xlab = 'Lesser hazard         Greater hazard                        ',
           vertices = TRUE,
           col=fpColors(box="deepskyblue4",line="deepskyblue4"))
}

pathway <- c("ABC_Transporter","Amino_Acid_Metabolism","Aminosugar","Cholesterol", "Cofactor","Complex_I","Complex_IV","Detox","Fatty_Acid","Glutathione","Glycan","Glycan_Anchor","Glycan_Degradation","Glycan_Sulfate","Glycolysis","Heparin_Sulfate","Hormone","Inositol_Phosphate","Ion_Transport","Krebs","Membrane_Lipid","Multipurpose","NAD","Neurotransmitter","Nucleotide", "Other","Other_Transport","Proton_Transport","Purine","Pyrimidine","Redox","Signalling","Small_Molecule)_Transport","Sphingolipid","Steroid","Sugar","Sulfate","Vitamin_A")

input <- scoresur
see(input)
input$pro[1:794]='ok'

data <- lapply(pathway,function(x)Subgroup_survival_analysis(pdata = input,
                                                          time ="os", status = "os_status",
                                                          variable = "pro",object = x))
foredata <- data[[1]][1,]
for (i in 2:38) {
  foredata <- rbind(foredata,data[[i]][1,])
}
rownames(foredata) <- colnames(input)[1:38]
foredata

rownames(foredata) <- c("ABC Transporter","Amino Acid Metabolism","Aminosugar","Cholesterol", "Cofactor","Complex I","Complex IV","Detox","Fatty Acid","Glutathione","Glycan","Glycan Anchor","Glycan Degradation","Glycan Sulfate","Glycolysis","Heparin Sulfate","Hormone","Inositol Phosphate","Ion Transport","Krebs","Membrane Lipid","Multipurpose","NAD","Neurotransmitter","Nucleotide", "Other","Other Transport","Proton Transport","Purine","Pyrimidine","Redox","Signalling","Small Molecule Transport","Sphingolipid","Steroid","Sugar","Sulfate","Vitamin A")
foredat <- as.matrix(foredata)
foredat <- cbind(as.data.frame(foredat),mean=apply(foredat[,3:4], 1,mean))
foredat <- foredat[order(foredat$mean),]
colnames(foredat)[3:4]=c('lower','upper')
cochrane_from_rmeta <- rbind(c(NA,NA,NA),foredat[,c(5,3,4)])

tabletext<-cbind(c("metabolic pathway",rownames(foredat)),
                 c('OR',foredat$HR),
                 c('95% CI',paste0('(',foredat$lower,',',foredat$upper,')'))
)

forestplot(tabletext, 
           graph.pos = 2,
           hrzl_lines = gpar(col="#444444"),
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,50),TRUE),
           clip=c(0.05,100), 
           xlog=TRUE,
           title = 'Prognostic associations of metabolic pathways',
           xlab = ' Lesser hazard         Greater hazard',
           vertices = TRUE,
           col=fpColors(box="deepskyblue4",line="deepskyblue4"))
