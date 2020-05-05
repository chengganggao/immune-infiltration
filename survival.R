rm(list = ls())  
library(survival)
library(survminer)

pan_sur <- pansurcluster
colnames(pan_sur)
pan_sur1 <-read.csv('pansur_0.05.csv')
pan_sur <- pan_sur1
pan_sur <- pan_sur[order(pan_sur$Monocytes),]
#pan_sur <- pan_sur[which(pan_sur$P.value<0.05),]
pan_sur$Monocytes
pan_sur$Monocytes_level <- NULL
pan_sur$Monocytes_level[127:482]='high'
pan_sur$Monocytes_level[1:126]='low'
pan_sur$Monocytes_level=as.factor(pan_sur$Monocytes_level)
pan_sur$Monocytes_level

fit <- survfit(Surv(os, os_status) ~ Monocytes_level, data = pan_sur)

ggsurvplot(
  fit, # survfit object with calculated statistics.
  pval = TRUE, # show p-value of log-rank test.
  conf.int = TRUE, # show confidence intervals for point estimaes of survival curves.
  conf.int.style = "ribbon", # customize style of confidence intervals
  xlab = "Time in months", # customize X axis label.
  break.time.by = 30, # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct", # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table
  ncensor.plot = TRUE, # plot the number of censored subjects at time t
  surv.median.line = "hv", # add the median survival pointer.
  legend.labs =c('Monocytes (high)','Monocytes (low)'), # change legend labels
  palette = c("#E7B800", "#2E9FDF")) # custom color palettes.




see(pdpdac)
pan_sur <- pansur
colnames(pan_sur)
pan_sur <- pan_sur[order(pan_sur$Monocytes),]
#pan_sur <- pan_sur[which(pan_sur$P.value<0.1),]
pan_sur$Monocytes
pan_sur$Monocytes_level <- NULL
pan_sur$Monocytes_level[428:794]='3'
pan_sur$Monocytes_level[191:427]='2'
pan_sur$Monocytes_level[1:190]='1'
pan_sur$Monocytes_level=as.factor(pan_sur$Monocytes_level)

pdsur <- pansurcluster[1:25,]
table(pdsur$nstage)
pan_sur <- rbind(pdsur[which(pdsur$nstage=='R0'),],pdsur[which(pdsur$nstage=='R1'),])
fit <- survfit(Surv(os, os_status) ~ nstage, data = pan_sur)

ggsurvplot(
  fit, # survfit object with calculated statistics.
  pval = TRUE, # show p-value of log-rank test.
  conf.int = TRUE, # show confidence intervals for point estimaes of survival curves.
  conf.int.style = "ribbon", # customize style of confidence intervals
  xlab = "Time in months", # customize X axis label.
  break.time.by = 30, # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct", # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table
  ncensor.plot = TRUE, # plot the number of censored subjects at time t
  surv.median.line = "hv", # add the median survival pointer.
  legend.labs =c('Negative resection margin','positive resection margin'), # change legend labels
  palette = c("#E7B800", "#2E9FDF")) # custom color palettes.


pan_sur_all_Monocytes <- pan_sur
pan_sur <- pan_sur_all_Monocytes



pan_sur <- pan_sur_all_Monocytes
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage=='male'),which(pan_sur$nstage=='female')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'G1'),which(pan_sur$nstage == 'G2'),which(pan_sur$nstage == 'G3')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'N0'),which(pan_sur$nstage == 'N1')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'T1'),which(pan_sur$nstage == 'T2'),which(pan_sur$nstage == 'T3'),which(pan_sur$nstage == 'T4')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'M0'),which(pan_sur$nstage == 'M1')),]
pan_sur_nstage$Monocytes_level

fit2 <- survfit( Surv(os, os_status) ~ Monocytes_level + nstage ,data = pan_sur_nstage)
ggsurvplot(fit2, conf.int = TRUE,
                     ggtheme = theme_bw())




sam1 <- sample(1:608,446,replace = F)
sam <- 1:893
sam[1:893] <- 'x'
sam[sam1]=NA
sam2 <- which(sam2=='x')
pan_sur$Monocytes_level[sam1]='low'
pan_sur$Monocytes_level[sam2]='high'
pan_sur$Monocytes_level=as.factor(pan_sur$Monocytes_level)


{ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE,
            risk.table = TRUE, # Add risk tablerisk.table.col = "strata", 
            # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            surv.median.line = "hv", # Specify median survival
            ggtheme = theme_bw(), # Change ggplot2 theme
            palette = c("#E7B800", "#2E9FDF"))}


pan_sur <- scoresur
colnames(pan_sur)
pan_sur <- pan_sur[order(pan_sur$Monocytes),]
#pan_sur <- pan_sur[which(pan_sur$P.value<0.1),]
pan_sur$Monocytes
pan_sur$Monocytes_level <- NULL
pan_sur$Monocytes_level[637:794]='3'
pan_sur$Monocytes_level[395:636]='2'
pan_sur$Monocytes_level[1:394]='1'
pan_sur$Monocytes_level=as.factor(pan_sur$Monocytes_level)



table(pdsur$nstage)
pan_sur <- rbind(pdsur[which(pdsur$nstage=='G1'),],pdsur[which(pdsur$nstage=='G2'),],pdsur[which(pdsur$nstage=='G3'),])

pansurcluster <- readcsv('pansurcluster.csv')
table(pansurcluster$clusterorder)
fit <- survfit(Surv(os, os_status) ~ clusterorder, data = pansurcluster)
see(pansurcluster)
write.csv(pansurcluster,'pansurcluster.csv')
identical(rownames(pansurcluster),rownames(clusterbar))
pansurcluster$clusterorder=clusterbar$clusterorder
ggsurvplot(
  fit, # survfit object with calculated statistics.
  pval = TRUE, # show p-value of log-rank test.
  conf.int = TRUE, # show confidence intervals for point estimaes of survival curves.
  conf.int.style = "ribbon", # customize style of confidence intervals
  xlab = "Time in months", # customize X axis label.
  break.time.by = 30, # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct", # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table
  ncensor.plot = TRUE, # plot the number of censored subjects at time t
  surv.median.line = "hv", # add the median survival pointer.
  legend.labs =c("cluster1","cluster2",'cluster3','cluster4','cluster5'), # change legend labels
  palette =c('#00DB00',"#9393FF","#FF5151","#E7B800",'yellow')) # custom color palettes.

"#E7B800"

score_Monocytes <- pan_sur




pan_sur <- pan_sur_all_Monocytes
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage=='male'),which(pan_sur$nstage=='female')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'G1'),which(pan_sur$nstage == 'G2'),which(pan_sur$nstage == 'G3')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'N0'),which(pan_sur$nstage == 'N1')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'T1'),which(pan_sur$nstage == 'T2'),which(pan_sur$nstage == 'T3'),which(pan_sur$nstage == 'T4')),]
pan_sur_nstage <- pan_sur[c(which(pan_sur$nstage == 'M0'),which(pan_sur$nstage == 'M1')),]
pan_sur_nstage$Monocytes_level

fit2 <- survfit( Surv(os, os_status) ~ Monocytes_level + nstage ,data = pan_sur_nstage)

dat <- as.data.frame(t(dat))
see(dat)

dat$os=as.numeric(dat$os)
dat$os_status[which(dat$os_status=='censored')]='0'
class(dat)




immgene <- immunegene[order(immunegene$Monocytes),]
gene <- immgene[c(1:500,17467:17966),]
gene <- pick(gene,dat)
rownames(gene)


R=immunegene['MCEMP1','Monocytes']
P=immuneP['THEMIS','Monocytes']
ymax=range(dat1$THEMIS)[2]
xmax=range(dat1$Monocytes)[2]

ggplot(data=dat1,aes(x=Monocytes,y=THEMIS))+geom_point()+geom_smooth()+geom_text(x = xmax-0.1, y = ymax, label = paste0('R=',R) ,parse = F)+geom_text(x = xmax-0.1, y = ymax-0.5, label = paste0('P=',P) ,parse = F)

dat <- readcsv('xtile3.csv')
dat$gene_level=NULL
dat <- dat[order(dat$HLA_B),]
dat$gene_level[179:235]='high'
dat$gene_level[1:178]='low'

fit <- survfit(Surv(os, os_status) ~ gene_level, data = dat)
ggsurvplot(fit, conf.int = TRUE,pval = TRUE,xlab = "Time in months",surv.median.line = "hv",
           legend.labs =c('HLA_B (high=57)','HLA_B (low=178)'),
           ggtheme = theme_bw())

