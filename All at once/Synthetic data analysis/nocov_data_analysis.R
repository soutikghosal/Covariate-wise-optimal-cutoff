

######################################
####
####  Analysis to perform analysis
####  without covariate
####
######################################

rm(list = ls())
library("rjags")
library(edgeR)
library(dplyr)
library(plyr)
library(ggplot2)
library(coda)
library(gtools)
library(reshape)
library(pracma)
library(mvtnorm)
library(readxl)

source("data_func.R")

#############
## Calling the data
##############

ADSP_data = readRDS("ADSPdata.rds")

######################
##  Calculating Table 1
######################

ADSP_data_sub = ADSP_data
myVars = colnames(ADSP_data_sub)
catVars = c("PHC_Sex","PHC_Race")

library(tableone)
ADSP_data_sub2 = ADSP_data_sub
myVars = colnames(ADSP_data_sub2)
catVars = c("PHC_Sex","PHC_Race")
tab.overall <- CreateTableOne(vars = myVars, data = ADSP_data_sub2, factorVars = catVars)
tab.Eye.cat <- CreateTableOne(vars = myVars, strata = "PHC_Diagnosis",
                              data = ADSP_data_sub2, factorVars = catVars)
tab.overall_p <- print(tab.overall, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tab.Eye.cat_p <- print(tab.Eye.cat, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tab.all = cbind(tab.overall_p,tab.Eye.cat_p)
tab.all = tab.all[-c(2:3),-c(3,ncol(tab.all))]
# write.csv(tab.all, file = "Table7.csv")

######################
##  Model fitting
######################

cutoff.summ.tab = function(tab, neg = TRUE){
  # tab = emp.fit$cutoff
  # if neg == TRUE, means biomarkers were multiplied by (-1)
  # if neg == FALSE, means biomarkers were raw, i.e. not multiplied by (-1)
  
  if(neg == TRUE){
    tab = tab[,c(1,2,4,3)]
    tab[,c(1,3,4)] = (-1)*tab[,c(1,3,4)]
  }
  out = paste0(round(tab[,1], 3), " (", 
               round(tab[,3], 3), ", ",
               round(tab[,4], 3), ")")
  return(out)
}

nocov_analysis = function(y0, y1, neg = TRUE){
  
  emp.fit = Empirical_Boot(y0, y1, B = 1000)
  NonPar.fit = NonPar_Boot(y0, y1, B = 1000)
  BN.fit = BiNormal_func(y0, y1, burnin = 5000, numb.iter = 5000)
  PV.fit = PV_func(y0, y1, burnin = 5000, numb.iter = 5000)
  SemiPV.fit = SemiPar_func(y0, y1, burnin = 5000, numb.iter = 5000)
  
  AUC.tab = rbind(emp.fit$AUC, NonPar.fit$AUC, BN.fit$AUC, 
                  PV.fit$AUC, SemiPV.fit$AUC)
  rownames(AUC.tab) = c("Empirical", "NonPar", "BN", "PV", "SemiPV")
  
  ROC.tab = cbind(emp.fit$ROC, NonPar.fit$ROC, BN.fit$ROC,
                  PV.fit$ROC, SemiPV.fit$ROC)
  colnames(ROC.tab) = c("Empirical", "NonPar", "BN", "PV", "SemiPV")
  
  cutoff.tab = rbind(cutoff.summ.tab(emp.fit$cutoff, neg),
                     cutoff.summ.tab(NonPar.fit$cutoff, neg),
                     cutoff.summ.tab(BN.fit$cutoff, neg),
                     cutoff.summ.tab(PV.fit$cutoff, neg),
                     cutoff.summ.tab(SemiPV.fit$cutoff, neg))
  colnames(cutoff.tab) = c("J","ER","CZ","IU")
  rownames(cutoff.tab) = c("Empirical", "NonPar", "BN", "PV", "SemiPV")
  
  out = list(AUC = AUC.tab, ROC = ROC.tab, cutoff = cutoff.tab)
  return(out)
}

##################
## AB42 analysis
##################

ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_AB42) %>% na.omit()
y0 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "3")]
summary(y0)
summary(y1)

y0 = -y0;
y1 = -y1;

tt = Sys.time()
AB42_nocov.fit = nocov_analysis(y0, y1, neg = TRUE)
Sys.time() - tt

grid = seq(0,1,length=100)

##################
## Tau analysis
##################

ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_Tau) %>% na.omit()
y0 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
summary(y0)
summary(y1)

tt = Sys.time()
tau_nocov.fit = nocov_analysis(y0, y1, neg = FALSE)
Sys.time() - tt


##################
## pTau analysis
##################

ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_pTau) %>% na.omit()
y0 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
summary(y0)
summary(y1)

tt = Sys.time()
ptau_nocov.fit = nocov_analysis(y0, y1, neg = FALSE)
Sys.time() - tt

# save.image("ADNI_ADSP_PHC_nocov.RData")
# load("ADNI_ADSP_PHC_nocov.RData")

######################
## Summarizing
######################

## AUC summary
AUC.summ.tab = function(AUC.tab){
  
  AUC.hat = paste0(round(AUC.tab[,1], 3),
                   " ",
                   paste0("(",round(AUC.tab[,3], 3), ", ",
                          round(AUC.tab[,4], 3),")"))
  tab = data.frame(Method = rownames(AUC.tab),
                   AUC = AUC.hat)
  colnames(tab)[2] = "AUC (95% CI)"
  
  return(tab)
}

## Cutoff summary
AUC.cutoff.summ.tab = function(fit){
  # AUC.tab = AB42_nocov.fit$AUC
  
  AUC.tab = fit$AUC
  AUC.hat = paste0(round(AUC.tab[,1], 3),
                   " ",
                   paste0("(",round(AUC.tab[,3], 3), ", ",
                          round(AUC.tab[,4], 3),")"))
  tab = data.frame(Method = rownames(AUC.tab),
                   AUC = AUC.hat)
  colnames(tab)[2] = "AUC (95% CI)"
  
  fin.AUC.tab = tab
  cutoff.tab = fit$cutoff
  
  fin.tab = data.frame(fin.AUC.tab, cutoff.tab)
  colnames(fin.tab)[2:6] = c("AUC (95% CI)", "J (95% CI)", "ER (95% CI)", "CZ (95% CI)",
                             "IU (95% CI)")
  return(fin.tab)
}

## ROC summary
ROC.summ.tab = function(ROC.tab){
  # ROC.tab = AB42_nocov.fit$ROC
  
  grid = seq(0,1, length = 100)
  ROC.tab = data.frame(grid,ROC.tab)
  ROC.melt = melt(ROC.tab, "grid")
  colnames(ROC.melt)[2] = c("Method") 
  
  cols = c("#999999", "#666666", "#000000",
           # "#006666", "#339999",
           "#CC0000", "#990000")
  dis.col = c("#66FF99","#FF9966")
  lty.type = c("twodash","longdash","dashed", #"1F","F1", 
               "4C88C488", "12345678")
  # c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  
  sz = 16
  plotROC = ggplot(ROC.melt, aes(x = grid, y = value, color = Method, group = Method)) +
    geom_line(aes(linetype=Method)) +
    xlab("1 - Specificity") + ylab("Sensitivity") +
    scale_color_manual(values = cols) +
    scale_linetype_manual(values = lty.type) +
    theme_classic() +
    theme(axis.text=element_text(size=sz),
          axis.title=element_text(size=sz),
          legend.position = "bottom",
          legend.title=element_text(size=sz),
          legend.text=element_text(size=sz))
  
  return(plotROC)
}

AUC_cutoff_tab = rbind(AUC.cutoff.summ.tab(AB42_nocov.fit),
                       AUC.cutoff.summ.tab(tau_nocov.fit),
                       AUC.cutoff.summ.tab(ptau_nocov.fit))
# write.csv(AUC_cutoff_tab, "Table8.csv")


## Plotting ROC together
require(ggpubr)
grid = seq(0,1,length=100)
all.ROC.plot = ggarrange(ROC.summ.tab(AB42_nocov.fit$ROC), 
                         ROC.summ.tab(tau_nocov.fit$ROC), 
                         ROC.summ.tab(ptau_nocov.fit$ROC), ncol = 3, nrow = 1,
                         labels = c("AB42","Tau","pTau"),
                         common.legend = TRUE, legend = "bottom")
# ggsave(filename = "Figure1.png",
#        plot = all.ROC.plot,
#        width = 12, height = 6,
#        device='png' #, dpi=1000
# )

###############
## Plotting biomarker densities and cutoffs
###############

## Function to summarize the plots
whole.summ = function(biom, y0,y1, fit){
  
  AUC.tab  = fit$AUC
  
  J.type = fit$cutoff[, 1]
  ER.type = fit$cutoff[, 2]
  CZ.type = fit$cutoff[, 3]
  IU.type = fit$cutoff[, 4]
  
  J = data.frame(value = sapply(J.type, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))), 
                 Method = rownames(AUC.tab))
  ER = data.frame(value = sapply(ER.type, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))), 
                  Method = rownames(AUC.tab))
  CZ = data.frame(value = sapply(CZ.type, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))), 
                  Method = rownames(AUC.tab))
  IU = data.frame(value = sapply(IU.type, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))), 
                  Method = rownames(AUC.tab))
  
  # J = data.frame(value = as.numeric(substr(fit$cutoff[,1],1,6)), Method = rownames(AUC.tab))
  # ER = data.frame(value = as.numeric(substr(fit$cutoff[,2],1,6)), Method = rownames(AUC.tab))
  # CZ = data.frame(value = as.numeric(substr(fit$cutoff[,3],1,6)), Method = rownames(AUC.tab))
  # IU = data.frame(value = as.numeric(substr(fit$cutoff[,4],1,6)), Method = rownames(AUC.tab))
  
  cols = c("#999999", "#666666", "#000000",
           # "#006666", "#339999",
           "#CC0000", "#990000")
  dis.col = c("#66FF99","#FF9966")
  
  sub.data = data.frame(biom = c(y0,y1), Category = c(rep("Normal cognition",length(y0)),
                                                      rep("AD",length(y1))))
  colnames(sub.data)[1] = biom
  sub.data$Category = factor(sub.data$Category, levels = c("Normal cognition","AD"))
  
  sz = 16
  
  J.plot = ggplot(sub.data,aes(x=sub.data[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = J,
      mapping = aes(xintercept = value, color = Method),
      linetype = "dashed", linewidth = 1, alpha = 0.99
    ) +
    scale_color_manual(values = cols)+
    scale_fill_manual(breaks = c("Normal cognition","AD"),values = dis.col)+
    ylab("Density") + xlab(biom) +
    theme_classic() +
    theme(axis.text=element_text(size=sz),
          axis.title=element_text(size=sz),
          legend.position = "bottom",
          legend.title=element_text(size=sz),
          legend.text=element_text(size=sz))
  ER.plot = ggplot(sub.data,aes(x=sub.data[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = ER,
      mapping = aes(xintercept = value, color = Method),
      linetype = "dashed", linewidth = 1, alpha = 0.99
    ) +
    scale_color_manual(values = cols)+
    scale_fill_manual(breaks = c("Normal cognition","AD"),values = dis.col)+
    ylab("Density") + xlab(biom) +
    theme_classic() +
    theme(axis.text=element_text(size=sz),
          axis.title=element_text(size=sz),
          legend.position = "bottom",
          legend.title=element_text(size=sz),
          legend.text=element_text(size=sz))
  CZ.plot = ggplot(sub.data,aes(x=sub.data[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = CZ,
      mapping = aes(xintercept = value, color = Method),
      linetype = "dashed", linewidth = 1, alpha = 0.99
    ) +
    scale_color_manual(values = cols)+
    scale_fill_manual(breaks = c("Normal cognition","AD"),values = dis.col)+
    ylab("Density") + xlab(biom) +
    theme_classic() +
    theme(axis.text=element_text(size=sz),
          axis.title=element_text(size=sz),
          legend.position = "bottom",
          legend.title=element_text(size=sz),
          legend.text=element_text(size=sz))
  IU.plot = ggplot(sub.data,aes(x=sub.data[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = IU,
      mapping = aes(xintercept = value, color = Method),
      linetype = "dashed", linewidth = 1, alpha = 0.99
    ) +
    scale_color_manual(values = cols)+
    scale_fill_manual(breaks = c("Normal cognition","AD"),values = dis.col)+
    ylab("Density") + xlab(biom) +
    theme_classic() +
    theme(axis.text=element_text(size=sz),
          axis.title=element_text(size=sz),
          legend.position = "bottom",
          legend.title=element_text(size=sz),
          legend.text=element_text(size=sz))
  
  require(ggpubr)
  all.cutoff.plot = ggarrange(J.plot, ER.plot, CZ.plot, IU.plot, ncol = 4, nrow = 1,
                              labels = c("J","ER","CZ","IU"),
                              common.legend = TRUE, legend = "bottom")
  out = list(J.plot= J.plot, ER.plot=ER.plot, CZ.plot=CZ.plot, IU.plot=IU.plot,
             all.cutoff.plot=all.cutoff.plot)
  return(out)
}

## PHC_AB42
ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_AB42) %>% na.omit()
y0 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "3")]
AB42.whole.summ = whole.summ(biom = "AB42", y0,y1, fit = AB42_nocov.fit)

## PHC_Tau
ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_Tau) %>% na.omit()
y0 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
Tau.whole.summ = whole.summ(biom = "Tau", y0,y1, fit = tau_nocov.fit)

## PHC_pTau
ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_pTau) %>% na.omit()
y0 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
pTau.whole.summ = whole.summ(biom = "pTau", y0,y1, fit = ptau_nocov.fit)

#### Cutoff plot saving

all.cutoff.plot = ggarrange(AB42.whole.summ$all.cutoff.plot,
                            Tau.whole.summ$all.cutoff.plot,
                            pTau.whole.summ$all.cutoff.plot,
                            ncol = 1, nrow = 3,
                            # labels = do.call(paste,expand.grid(bioms, sexs))[c(1,4,2,5,3,6)],
                            common.legend = TRUE, legend = "bottom")
# ggsave(filename = "Figure2.png",
#        plot = all.cutoff.plot,
#        width = 18, height = 12,
#        device='png' #, dpi=1000
# )
