
######################################
####
####  Analysis to perform analysis
####  with covariate "Sex".
####  Sex-specific estimates
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

source("data_func_cov.R")

#############
## Calling the data
##############

ADSP_data = readRDS("ADSPdata.rds")

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

cov_analysis = function(y0, y1, X0, X1, x = c(0,1), neg = TRUE){
  
  BN.fit = BiNormal_reg_fit(y0, y1, X0, X1, x = c(0,1), burnin = 5000, numb.iter = 5000)
  PV.fit = BayesPV_reg_fit(y0, y1, X0, X1, x = c(0,1), burnin = 5000, numb.iter = 5000)
  SemiPV.fit = SemiPar_reg_fit(y0, y1, X0, X1, x = c(0,1), burnin = 5000, numb.iter = 5000)
  
  AUC0.tab = rbind(BN.fit$AUC$AUC0, PV.fit$AUC$AUC0, SemiPV.fit$AUC$AUC0)
  AUC1.tab = rbind(BN.fit$AUC$AUC1, PV.fit$AUC$AUC1, SemiPV.fit$AUC$AUC1)
  rownames(AUC0.tab) = rownames(AUC1.tab) = c("BN", "PV", "SemiPV")
  
  ROC0.tab = cbind(BN.fit$ROC$ROC0, PV.fit$ROC$ROC0, SemiPV.fit$ROC$ROC0)
  ROC1.tab = cbind(BN.fit$ROC$ROC1, PV.fit$ROC$ROC1, SemiPV.fit$ROC$ROC1)
  colnames(ROC0.tab) = colnames(ROC1.tab) = c("BN", "PV", "SemiPV")
  
  cutoff0.tab = rbind(cutoff.summ.tab(BN.fit$cutoff$cutoff0, neg),
                      cutoff.summ.tab(PV.fit$cutoff$cutoff0, neg),
                      cutoff.summ.tab(SemiPV.fit$cutoff$cutoff0, neg))
  cutoff1.tab = rbind(cutoff.summ.tab(BN.fit$cutoff$cutoff1, neg),
                      cutoff.summ.tab(PV.fit$cutoff$cutoff1, neg),
                      cutoff.summ.tab(SemiPV.fit$cutoff$cutoff1, neg))
  colnames(cutoff0.tab) = colnames(cutoff1.tab) = c("J","ER","CZ","IU")
  rownames(cutoff0.tab) = rownames(cutoff1.tab) = c("BN", "PV", "SemiPV")
  
  out = list(AUC0 = AUC0.tab, ROC0 = ROC0.tab, cutoff0 = cutoff0.tab,
             AUC1 = AUC1.tab, ROC1 = ROC1.tab, cutoff1 = cutoff1.tab)
  return(out)
}

check_data = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_AB42, PHC_Sex, PHC_Race, PHC_AGE) %>%
  na.omit()
check_data0 = check_data %>% dplyr::filter(PHC_Diagnosis == "1")
check_data1 = check_data %>% dplyr::filter(PHC_Diagnosis == "3")
X0 = check_data0 %>% dplyr::select(PHC_Sex, PHC_Race, PHC_AGE)
X1 = check_data1 %>% dplyr::select(PHC_Sex, PHC_Race, PHC_AGE)

### Modifying the SEX variable

ADSP_data$PHC_Sex_New = ADSP_data$PHC_Sex
ADSP_data$PHC_Sex_New[which(ADSP_data$PHC_Sex == "1")] = 0  ## Male
ADSP_data$PHC_Sex_New[which(ADSP_data$PHC_Sex == "2")] = 1  ## Female
ADSP_data$PHC_Sex_New = as.numeric(ADSP_data$PHC_Sex_New)

##################
## AB42 analysis
##################

ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_AB42, PHC_Sex_New) %>%
  na.omit()

y0 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "3")]
X0 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "1")]
X1 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "3")]
y0 = -y0;
y1 = -y1

tt = Sys.time()
AB42_cov.fit = cov_analysis(y0, y1, X0, X1, neg = TRUE)
Sys.time() - tt


##################
## Tau analysis
##################

ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_Tau, PHC_Sex_New) %>% na.omit()
y0 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
X0 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "1")]
X1 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "3")]
summary(y0)
summary(y1)

tt = Sys.time()
tau_cov.fit = cov_analysis(y0, y1, X0, X1, neg = FALSE)
Sys.time() - tt


##################
## pTau analysis
##################

ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_pTau, PHC_Sex_New) %>% na.omit()
y0 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
X0 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "1")]
X1 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "3")]
summary(y0)
summary(y1)

tt = Sys.time()
ptau_cov.fit = cov_analysis(y0, y1, X0, X1, neg = FALSE)
Sys.time() - tt


# save.image("ADNI_ADSP_PHC_sex.RData")
# load("ADNI_ADSP_PHC_sex.RData")


######################
## Summarizing
######################

AUC.summ.tab = function(AUC0, AUC1){
  # AUC.tab = AB42_nocov.fit$AUC
  
  AUC.hat0 = paste0(round(AUC0[,1], 3),
                    " ",
                    paste0("(",round(AUC0[,3], 3), ", ",
                           round(AUC0[,4], 3),")"))
  AUC.hat1 = paste0(round(AUC1[,1], 3),
                    " ",
                    paste0("(",round(AUC1[,3], 3), ", ",
                           round(AUC1[,4], 3),")"))
  
  tab = data.frame(Method = rownames(AUC0),
                   AUC0 = AUC.hat0,
                   AUC1 = AUC.hat1)
  colnames(tab)[2:3] = c("AUC0 (95% CI)","AUC1 (95% CI)")
  return(tab)
}

AUC.cutoff.summ.tab = function(fit){
  # AUC.tab = AB42_cov.fit$AUC
  
  AUC0 = fit$AUC0
  AUC1 = fit$AUC1
  tab = AUC.summ.tab(AUC0, AUC1)
  
  cutoff.tab0 = fit$cutoff0
  cutoff.tab1 = fit$cutoff1
  
  fin.tab0 = data.frame(tab[,-3], cutoff.tab0)
  fin.tab1 = data.frame(tab[,-2], cutoff.tab1)
  rownames(fin.tab0) = rownames(fin.tab1) = NULL
  colnames(fin.tab0)[2:6] = colnames(fin.tab1)[2:6] = c("AUC (95% CI)", "J (95% CI)", "ER (95% CI)", "CZ (95% CI)",
                                                        "IU (95% CI)")
  out = list(fin.tab0 = fin.tab0, fin.tab1 = fin.tab1)
  return(out)
}

ROC.summ.tab = function(ROC0, ROC1){
  
  grid = seq(0,1, length = 100)
  ROC.tab0 = data.frame(grid,ROC0)
  ROC.tab1 = data.frame(grid,ROC1)
  
  ROC.melt0 = melt(ROC.tab0, "grid")
  ROC.melt1 = melt(ROC.tab1, "grid")
  colnames(ROC.melt0)[2] = colnames(ROC.melt1)[2] = c("Method") 
  
  cols = c("#000000", "#CC0000", "#990000")
  dis.col = c("#66FF99","#FF9966")
  lty.type = c("dashed", "4C88C488", "12345678")
  # c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  
  sz = 16
  
  plotROC0 = ggplot(ROC.melt0, aes(x = grid, y = value, color = Method, group = Method)) +
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
  
  plotROC1 = ggplot(ROC.melt1, aes(x = grid, y = value, color = Method, group = Method)) +
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
  
  out = list(plotROC0=plotROC0, plotROC1=plotROC1)
  return(out)
}


AB42.cutoff.tab = rbind(AUC.cutoff.summ.tab(AB42_cov.fit)$fin.tab0,
                        AUC.cutoff.summ.tab(AB42_cov.fit)$fin.tab1)
Tau.cutoff.tab = rbind(AUC.cutoff.summ.tab(tau_cov.fit)$fin.tab0,
                       AUC.cutoff.summ.tab(tau_cov.fit)$fin.tab1)
pTau.cutoff.tab = rbind(AUC.cutoff.summ.tab(ptau_cov.fit)$fin.tab0,
                        AUC.cutoff.summ.tab(ptau_cov.fit)$fin.tab1)
AUC_cutoff_tab = rbind(AB42.cutoff.tab,Tau.cutoff.tab,pTau.cutoff.tab)
# write.csv(AUC_cutoff_tab, "Table9.csv")


## Plotting ROC together

require(ggpubr)
bioms = c("AB42","Tau","pTau")
sexs = c("(Sex: Male)", "(Sex: Female)")
all.ROC.plot = ggarrange(ROC.summ.tab(AB42_cov.fit$ROC0, AB42_cov.fit$ROC1)$plotROC0,
                         ROC.summ.tab(AB42_cov.fit$ROC0, AB42_cov.fit$ROC1)$plotROC1,
                         ROC.summ.tab(tau_cov.fit$ROC0, tau_cov.fit$ROC1)$plotROC0,
                         ROC.summ.tab(tau_cov.fit$ROC0, tau_cov.fit$ROC1)$plotROC1,
                         ROC.summ.tab(ptau_cov.fit$ROC0, ptau_cov.fit$ROC1)$plotROC0,
                         ROC.summ.tab(ptau_cov.fit$ROC0, ptau_cov.fit$ROC1)$plotROC1,
                         ncol = 2, nrow = 3,
                         labels = do.call(paste,expand.grid(bioms, sexs))[c(1,4,2,5,3,6)],
                         common.legend = TRUE, legend = "bottom")

# ggsave(filename = "Figure3.png",
#        plot = all.ROC.plot,
#        width = 12, height = 18,
#        device='png' #, dpi=1000
# )


###############
## Plotting biomarker densities and cutoffs
###############

## Function to summarize the plots
whole.summ = function(biom, y0,y1,X0,X1, fit){
  
  # fit = AB42_cov.fit
  
  AUC0 = fit$AUC0
  AUC1 = fit$AUC1
  tab = AUC.summ.tab(AUC0, AUC1)
  
  J.type0 = fit$cutoff0[, 1]
  ER.type0 = fit$cutoff0[, 2]
  CZ.type0 = fit$cutoff0[, 3]
  IU.type0 = fit$cutoff0[, 4]
  
  J.type1 = fit$cutoff1[, 1]
  ER.type1 = fit$cutoff1[, 2]
  CZ.type1 = fit$cutoff1[, 3]
  IU.type1 = fit$cutoff1[, 4]
  
  
  J = data.frame(value0 = sapply(J.type0, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                 value1 = sapply(J.type1, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                 Method = rownames(AUC0))
  ER = data.frame(value0 = sapply(ER.type0, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                  value1 = sapply(ER.type1, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                  Method = rownames(AUC0))
  CZ = data.frame(value0 = sapply(CZ.type0, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                  value1 = sapply(CZ.type1, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                  Method = rownames(AUC0))
  IU = data.frame(value0 = sapply(IU.type0, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                  value1 = sapply(IU.type1, function(x) as.numeric(sapply(strsplit(x, " "), "[[", 1))),
                  Method = rownames(AUC0))
  
  cols = c("#000000","#CC0000", "#006666")
  dis.col = c("#66FF99","#FF9966")
  
  sub.data = data.frame(biom = c(y0,y1), Sex = c(X0,X1), Category = c(rep("Normal cognition",length(y0)),
                                                                      rep("AD",length(y1))))
  colnames(sub.data)[1] = biom
  sub.data$Category = factor(sub.data$Category, levels = c("Normal cognition","AD"))
  
  sub.data0 = sub.data %>% filter(Sex==0)
  sub.data1 = sub.data %>% filter(Sex==1)
  
  sz = 16
  
  J.plot0 = ggplot(sub.data0,aes(x=sub.data0[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = J,
      mapping = aes(xintercept = value0, color = Method),
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
  
  J.plot1 = ggplot(sub.data1,aes(x=sub.data1[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = J,
      mapping = aes(xintercept = value1, color = Method),
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
  
  
  ER.plot0 = ggplot(sub.data0,aes(x=sub.data0[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = ER,
      mapping = aes(xintercept = value0, color = Method),
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
  ER.plot1 = ggplot(sub.data1,aes(x=sub.data1[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = ER,
      mapping = aes(xintercept = value1, color = Method),
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
  
  CZ.plot0 = ggplot(sub.data0,aes(x=sub.data0[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = CZ,
      mapping = aes(xintercept = value0, color = Method),
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
  CZ.plot1 = ggplot(sub.data1,aes(x=sub.data1[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = CZ,
      mapping = aes(xintercept = value1, color = Method),
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
  
  IU.plot0 = ggplot(sub.data0,aes(x=sub.data0[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = IU,
      mapping = aes(xintercept = value0, color = Method),
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
  IU.plot1 = ggplot(sub.data1,aes(x=sub.data1[,1], fill=Category)) + geom_density(alpha=0.25) +
    geom_vline(
      data = IU,
      mapping = aes(xintercept = value1, color = Method),
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
  all.cutoff.plot = ggarrange(J.plot0,J.plot1, ER.plot0,ER.plot1, 
                              CZ.plot0,CZ.plot1, IU.plot0,IU.plot1,
                              ncol = 2, nrow = 4,
                              labels = do.call(paste,expand.grid(c("J","ER","CZ","IU"), c("(Sex: Male)","(Sex: Female)")))[c(1,5,2,6,3,7,4,8)],
                              common.legend = TRUE, legend = "bottom")
  out = list(J.plot0= J.plot0, J.plot1=J.plot1, ER.plot0=ER.plot0, 
             ER.plot1=ER.plot1, CZ.plot0=CZ.plot0, CZ.plot1=CZ.plot1, 
             IU.plot0=IU.plot0, IU.plot1=IU.plot1,
             all.cutoff.plot=all.cutoff.plot)
  return(out)
}


## PHC_AB42
ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_AB42, PHC_Sex_New) %>%
  na.omit()
y0 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_AB42[which(ADSP_data_sub$PHC_Diagnosis == "3")]
X0 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "1")]
X1 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "3")]
AB42.whole.summ = whole.summ(biom = "AB42", y0,y1,X0,X1, fit = AB42_cov.fit)

## PHC_Tau
ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_Tau, PHC_Sex_New) %>% na.omit()
y0 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_Tau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
X0 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "1")]
X1 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "3")]
Tau.whole.summ = whole.summ(biom = "Tau", y0,y1,X0,X1, fit = tau_cov.fit)

# PHC_pTau
ADSP_data_sub = ADSP_data %>% dplyr::select(PHC_Diagnosis, PHC_pTau, PHC_Sex_New) %>% na.omit()
y0 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "1")]
y1 = ADSP_data_sub$PHC_pTau[which(ADSP_data_sub$PHC_Diagnosis == "3")]
X0 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "1")]
X1 = ADSP_data_sub$PHC_Sex_New[which(ADSP_data_sub$PHC_Diagnosis == "3")]
pTau.whole.summ = whole.summ(biom = "pTau", y0,y1,X0,X1, fit = ptau_cov.fit)


#### Cutoff plot saving

# ggsave(filename = "Figure4.png",
#        plot = AB42.whole.summ$all.cutoff.plot,
#        width = 12, height = 18,
#        device='png' #, dpi=1000
# )
# 
# ggsave(filename = "Figure5.png",
#        plot = Tau.whole.summ$all.cutoff.plot,
#        width = 12, height = 18,
#        device='png' #, dpi=1000
# )
# 
# ggsave(filename = "Figure6.png",
#        plot = pTau.whole.summ$all.cutoff.plot,
#        width = 12, height = 18,
#        device='png' #, dpi=1000
# )

