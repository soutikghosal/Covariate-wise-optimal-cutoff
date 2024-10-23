

rm(list = ls())
# setwd("C:/Users/gph7ps/OneDrive - University of Virginia/Methodological Research/Optimal Cutoff ROC/Simulation")
setwd("~/Library/CloudStorage/OneDrive-UniversityofVirginia/Methodological Research/Optimal Cutoff ROC/Simulation")

library(reshape)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(ggpubr)


summary.tab = function(AUC, AUC.True, ROC, ROC.True, 
                       J, ER, CZ, IU, Cutoff.tab){
  
  AUC.hat = apply(AUC,2,mean.NA)
  AUC.sd = apply(AUC,2,sd.NA)
  AUC.bias = apply(AUC-AUC.True, 2,mean.NA)
  AUC.bias.SD = apply(AUC-AUC.True, 2,sd.NA)
  AUC.bias.whole = AUC-AUC.True
  AUC.mse = apply((AUC-AUC.True)^2, 2,mean.NA)
  EMSE = apply(apply((ROC-ROC.True)^2,c(2,3),mean.NA),2,mean.NA)
  AUC.tab = data.frame(True.AUC = AUC.True,
                       AUC = AUC.hat, SD = AUC.sd, Bias = AUC.bias, Bias.SD = AUC.bias.SD,
                       MSE = AUC.mse, EMSE = EMSE)
  rownames(AUC.tab) = c("Emp","BN","BG","BiChi","NonPar",
                        "PV", "Semi.PV")
  ROC.hat = apply(ROC,c(2,3),mean.NA)
  colnames(ROC.hat) = c("Emp","BN","BG","BiChi","NonPar",
                        "PV", "Semi.PV")
  
  J.hat = apply(J, c(2,3), mean.NA)
  ER.hat = apply(ER, c(2,3), mean.NA)
  CZ.hat = apply(CZ, c(2,3), mean.NA)
  IU.hat = apply(IU, c(2,3), mean.NA)
  
  Bias.J.hat = apply(J - Cutoff.tab[1,], c(2,3), mean.NA)
  Bias.ER.hat = apply(ER - Cutoff.tab[2,], c(2,3), mean.NA)
  Bias.CZ.hat = apply(CZ - Cutoff.tab[3,], c(2,3), mean.NA)
  Bias.IU.hat = apply(IU - Cutoff.tab[4,], c(2,3), mean.NA)
  
  Bias.J.SD = apply(J - Cutoff.tab[1,], c(2,3), sd.NA)
  Bias.ER.SD = apply(ER - Cutoff.tab[2,], c(2,3), sd.NA)
  Bias.CZ.SD = apply(CZ - Cutoff.tab[3,], c(2,3), sd.NA)
  Bias.IU.SD = apply(IU - Cutoff.tab[4,], c(2,3), sd.NA)
  
  J.cutoff.bias = J[,1,] - Cutoff.tab[1,1]
  ER.cutoff.bias = ER[,1,] - Cutoff.tab[2,1]
  CZ.cutoff.bias = CZ[,1,] - Cutoff.tab[3,1]
  IU.cutoff.bias = IU[,1,] - Cutoff.tab[4,1]
  
  colnames(AUC.bias.whole) = colnames(J.cutoff.bias) = colnames(ER.cutoff.bias) = 
    colnames(CZ.cutoff.bias) = colnames(IU.cutoff.bias) = c("Emp","BN","BG","BiChi","NonPar",
                                                            "PV", "Semi.PV")
  
  # apply(J.cutoff.bias,2,median.NA)
  # apply(J - Cutoff.tab[1,], c(2,3), median.NA)
  
  AUC.bias.2 = melt(AUC.bias.whole)
  J.cutoff.bias.2 = melt(J.cutoff.bias)
  ER.cutoff.bias.2 = melt(ER.cutoff.bias)
  CZ.cutoff.bias.2 = melt(CZ.cutoff.bias)
  IU.cutoff.bias.2 = melt(IU.cutoff.bias)
  AUC.bias.2$Method = "AUC"
  J.cutoff.bias.2$Method = "J"
  ER.cutoff.bias.2$Method = "ER"
  CZ.cutoff.bias.2$Method = "CZ"
  IU.cutoff.bias.2$Method = "IU"
  bias.whole = rbind(AUC.bias.2, J.cutoff.bias.2, ER.cutoff.bias.2,
                     CZ.cutoff.bias.2, IU.cutoff.bias.2)
  colnames(bias.whole)[2:3] = c("Model", "Bias")
  bias.whole$X1 = NULL
  
  colnames(J.hat) = colnames(ER.hat) = colnames(CZ.hat) = 
    colnames(IU.hat) = colnames(Bias.J.hat) = colnames(Bias.ER.hat) = colnames(Bias.CZ.hat) = 
    colnames(Bias.IU.hat) = c("Emp","BN","BG","BiChi","NonPar",
                              "PV", "Semi.PV")
  rownames(J.hat) = rownames(ER.hat) = rownames(CZ.hat) = 
    rownames(IU.hat) = rownames(Bias.J.hat) = rownames(Bias.ER.hat) = rownames(Bias.CZ.hat) = 
    rownames(Bias.IU.hat) = c("cutoff","value","sensitivity", "specificity", 
                              "sensitivity (data)", "specificity (data)")
  
  
  line_size <- 1 # defining variable upfront as we will re-use it
  base_size <- 16 # defining separately, same as for line_size
  axis_text_rel_size = -1
  title_text_rel_size = +2
  
  
  AUC.plot = ggplot(bias.whole %>% filter(Method == "AUC"), aes(x = Model, y=Bias)) +
    geom_boxplot()
  J.plot = ggplot(bias.whole %>% filter(Method == "J"), aes(x = Model, y=Bias)) +
    geom_boxplot()
  ER.plot = ggplot(bias.whole %>% filter(Method == "ER"), aes(x = Model, y=Bias)) +
    geom_boxplot()
  CZ.plot = ggplot(bias.whole %>% filter(Method == "CZ"), aes(x = Model, y=Bias)) +
    geom_boxplot()
  IU.plot = ggplot(bias.whole %>% filter(Method == "IU"), aes(x = Model, y=Bias)) +
    geom_boxplot()
  
  whole.plot = ggplot(bias.whole, aes(x = Model, y=Bias, fill = Method)) +
    geom_boxplot() + # ylim(c(-1, 1)) +
    geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) +
    scale_fill_brewer(palette="RdBu") +
    theme_foundation(base_size = base_size, base_family = "sans") + 
    theme(
      legend.title=element_blank(),
      legend.position = "bottom",
      # aspect.ratio = 1.6,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", 
                                size = rel((title_text_rel_size + base_size) / base_size), hjust = 0.5),
      axis.line = element_line(colour="black", linewidth = line_size),
      axis.ticks = element_line(colour="black", linewidth = line_size),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
      # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.background = element_blank())
  
  
  allbutAUC.plot = ggplot(bias.whole %>% filter(Method != "AUC"), aes(x = Model, y=Bias, fill = Method)) +
    geom_boxplot() + # ylim(c(-1, 1)) +
    geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) +
    scale_fill_brewer(palette="RdBu") +
    theme_foundation(base_size = base_size, base_family = "sans") + 
    theme(
      legend.title=element_blank(),
      legend.position = "bottom",
      # aspect.ratio = 1.6,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", 
                                size = rel((title_text_rel_size + base_size) / base_size), hjust = 0.5),
      axis.line = element_line(colour="black", linewidth = line_size),
      axis.ticks = element_line(colour="black", linewidth = line_size),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
      # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.background = element_blank())
  
  bias.whole.noCon = bias.whole %>% filter(Method != "AUC" & Model != "BG" & Model != "BiChi")
  bias.whole.noCon = droplevels(bias.whole.noCon)
  allbutAUC.noCon.plot = ggplot(bias.whole.noCon, aes(x = Model, y=Bias, fill = Method)) +
    geom_boxplot() + # ylim(c(-1, 1)) +
    geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) +
    scale_fill_brewer(palette="RdBu") +
    theme_foundation(base_size = base_size, base_family = "sans") + 
    theme(
      legend.title=element_blank(),
      legend.position = "bottom",
      # aspect.ratio = 1.6,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", 
                                size = rel((title_text_rel_size + base_size) / base_size), hjust = 0.5),
      axis.line = element_line(colour="black", linewidth = line_size),
      axis.ticks = element_line(colour="black", linewidth = line_size),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
      # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.background = element_blank())
  
  bias.noCon = bias.whole %>% filter(Model != "BG" & Model != "BiChi")
  AUC.noCon.plot = ggplot(bias.noCon %>% filter(Method == "AUC"), aes(x = Model, y=Bias)) +
    geom_boxplot() + # ylim(c(-1, 1)) +
    geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) +
    ylab("AUC") +
    theme_foundation(base_size = base_size, base_family = "sans") + 
    theme(
      legend.title=element_blank(),
      legend.position = "bottom",
      # aspect.ratio = 1.6,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", 
                                size = rel((title_text_rel_size + base_size) / base_size), hjust = 0.5),
      axis.line = element_line(colour="black", linewidth = line_size),
      axis.ticks = element_line(colour="black", linewidth = line_size),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
      # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.background = element_blank())
  
  
  out = list(AUC.tab = AUC.tab, ROC.hat = ROC.hat, True.ROC = ROC.True,
             J = J.hat, ER = ER.hat, CZ = CZ.hat,
             IU = IU.hat,
             Bias.J = Bias.J.hat, Bias.ER = Bias.ER.hat, Bias.CZ = Bias.CZ.hat,
             Bias.IU = Bias.IU.hat,
             Bias.J.SD = Bias.J.SD, Bias.ER.SD = Bias.ER.SD, Bias.CZ.SD = Bias.CZ.SD,
             Bias.IU.SD = Bias.IU.SD,
             AUC.bias.whole=AUC.bias.whole, J.cutoff.bias=J.cutoff.bias, 
             ER.cutoff.bias=ER.cutoff.bias, CZ.cutoff.bias=CZ.cutoff.bias, 
             IU.cutoff.bias=IU.cutoff.bias,
             bias.whole = bias.whole,
             AUC.plot=AUC.plot, J.plot=J.plot, ER.plot=ER.plot, CZ.plot=CZ.plot, 
             IU.plot=IU.plot, whole.plot = whole.plot, allbutAUC.plot = allbutAUC.plot,
             allbutAUC.noCon.plot = allbutAUC.noCon.plot,
             AUC.noCon.plot = AUC.noCon.plot)
  return(out)
}


bias.table = function(summ.object){
  
  library(dplyr)

  AUC.bias = list(Bias = apply(summ.object$AUC.bias.whole,2,mean.NA),
                  Bias.SD = apply(summ.object$AUC.bias.whole,2,sd.NA))
  
  # AUC.bias = summ.object$AUC.tab %>% select(Bias, Bias.SD)
  nms = names(AUC.bias$Bias)
  AUC.bias = paste(round(AUC.bias$Bias,3), "±", round(AUC.bias$Bias.SD,3), sep = "")
  
  # J.Bias = data.frame(Bias = summ.object$Bias.J[1,], Bias.SD = summ.object$Bias.J.SD[1,])
  J.Bias = data.frame(Bias = apply(summ.object$J.cutoff.bias,2,mean.NA), Bias.SD = apply(summ.object$J.cutoff.bias,2,sd.NA))
  J.Bias = paste(round(J.Bias$Bias,3), "±", round(J.Bias$Bias.SD,3), sep = "")
  
  # ER.Bias = data.frame(Bias = summ.object$Bias.ER[1,], Bias.SD = summ.object$Bias.ER.SD[1,])
  ER.Bias = data.frame(Bias = apply(summ.object$ER.cutoff.bias,2,mean.NA), Bias.SD = apply(summ.object$ER.cutoff.bias,2,sd.NA))
  ER.Bias = paste(round(ER.Bias$Bias,3), "±", round(ER.Bias$Bias.SD,3), sep = "")
  
  # CZ.Bias = data.frame(Bias = summ.object$Bias.CZ[1,], Bias.SD = summ.object$Bias.CZ.SD[1,])
  CZ.Bias = data.frame(Bias = apply(summ.object$CZ.cutoff.bias,2,mean.NA), Bias.SD = apply(summ.object$CZ.cutoff.bias,2,sd.NA))
  CZ.Bias = paste(round(CZ.Bias$Bias,3), "±", round(CZ.Bias$Bias.SD,3), sep = "")
  
  # IU.Bias = data.frame(Bias = summ.object$Bias.IU[1,], Bias.SD = summ.object$Bias.IU.SD[1,])
  IU.Bias = data.frame(Bias = apply(summ.object$IU.cutoff.bias,2,mean.NA), Bias.SD = apply(summ.object$IU.cutoff.bias,2,sd.NA))
  IU.Bias = paste(round(IU.Bias$Bias,3), "±", round(IU.Bias$Bias.SD,3), sep = "")
  
  final.bias = data.frame(AUC = AUC.bias, J = J.Bias, ER = ER.Bias,
                          CZ = CZ.Bias, IU = IU.Bias)
  rownames(final.bias) = nms
  return(final.bias)
}

bias.table.me = function(summ.object){
  
  median.NA = function(x){as.numeric(quantile(x, probs = c(0.5), na.rm = TRUE))}
  IQR.NA = function(x){IQR(x, na.rm = TRUE, type = 7)}
  
  library(dplyr)

  AUC.bias = list(Bias = apply(summ.object$AUC.bias.whole,2,median.NA),
                  Bias.IQR = apply(summ.object$AUC.bias.whole,2,IQR.NA))
  
  # AUC.bias = summ.object$AUC.tab %>% select(Bias, Bias.IQR)
  nms = names(AUC.bias$Bias)
  AUC.bias = paste(round(AUC.bias$Bias,3), "±", round(AUC.bias$Bias.IQR,3), sep = "")
  
  # J.Bias = data.frame(Bias = summ.object$Bias.J[1,], Bias.IQR = summ.object$Bias.J.IQR[1,])
  J.Bias = data.frame(Bias = apply(summ.object$J.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$J.cutoff.bias,2,IQR.NA))
  J.Bias = paste(round(J.Bias$Bias,3), "±", round(J.Bias$Bias.IQR,3), sep = "")
  
  # ER.Bias = data.frame(Bias = summ.object$Bias.ER[1,], Bias.IQR = summ.object$Bias.ER.IQR[1,])
  ER.Bias = data.frame(Bias = apply(summ.object$ER.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$ER.cutoff.bias,2,IQR.NA))
  ER.Bias = paste(round(ER.Bias$Bias,3), "±", round(ER.Bias$Bias.IQR,3), sep = "")
  
  # CZ.Bias = data.frame(Bias = summ.object$Bias.CZ[1,], Bias.IQR = summ.object$Bias.CZ.IQR[1,])
  CZ.Bias = data.frame(Bias = apply(summ.object$CZ.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$CZ.cutoff.bias,2,IQR.NA))
  CZ.Bias = paste(round(CZ.Bias$Bias,3), "±", round(CZ.Bias$Bias.IQR,3), sep = "")
  
  # IU.Bias = data.frame(Bias = summ.object$Bias.IU[1,], Bias.IQR = summ.object$Bias.IU.IQR[1,])
  IU.Bias = data.frame(Bias = apply(summ.object$IU.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$IU.cutoff.bias,2,IQR.NA))
  IU.Bias = paste(round(IU.Bias$Bias,3), "±", round(IU.Bias$Bias.IQR,3), sep = "")
  
  final.bias = data.frame(AUC = AUC.bias, J = J.Bias, ER = ER.Bias,
                          CZ = CZ.Bias, IU = IU.Bias)
  rownames(final.bias) = nms
  return(final.bias)
}


##################################
##################################
####                          ####
#### Low Sample size, Low AUC ####
####                          ####
##################################
##################################

Samp_Cat = "Low"
AUC_Cat = "Low"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")

dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_LowSamp_LowAUC.csv")


require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_lowSamp_lowAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_lowSamp_lowAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

#####################################
#####################################
####                             ####
#### Low Sample size, Medium AUC ####
####                             ####
#####################################
#####################################

Samp_Cat = "Low"
AUC_Cat = "Medium"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_LowSamp_MedAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])


ggsave(filename = "Plot_sim_bias_noCon_lowSamp_medAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_lowSamp_medAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

#####################################
#####################################
####                             ####
#### Low Sample size, High AUC   ####
####                             ####
#####################################
#####################################

Samp_Cat = "Low"
AUC_Cat = "High"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_", k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_", k, ".rds")
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_LowSamp_HighAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_lowSamp_highAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_lowSamp_highAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

#####################################
#####################################
####                             ####
#### Medium Sample size, Low AUC ####
####                             ####
#####################################
#####################################

Samp_Cat = "Medium"
AUC_Cat = "Low"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 

}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_MedSamp_LowAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_medSamp_lowAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_medSamp_lowAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

########################################
########################################
####                                ####
#### Medium Sample size, Medium AUC ####
####                                ####
########################################
########################################

Samp_Cat = "Medium"
AUC_Cat = "Medium"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",1, ".rds")

dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_MedSamp_MedAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])


ggsave(filename = "Plot_sim_bias_noCon_medSamp_medAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_medSamp_medAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

########################################
########################################
####                                ####
#### Medium Sample size, High AUC   ####
####                                ####
########################################
########################################

Samp_Cat = "Medium"
AUC_Cat = "High"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",1, ".rds")

dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_MedSamp_HighAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_medSamp_highAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_medSamp_highAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

###################################
###################################
####                           ####
#### High Sample size, Low AUC ####
####                           ####
###################################
###################################

Samp_Cat = "High"
AUC_Cat = "Low"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")

  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_HighSamp_LowAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-0.5, 0.5)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-1, 1)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-1,1)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_highSamp_lowAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_highSamp_lowAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

######################################
######################################
####                              ####
#### High Sample size, Medium AUC ####
####                              ####
######################################
######################################

Samp_Cat = "High"
AUC_Cat = "Medium"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)

bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_HighSamp_MedAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-0.5, 0.5)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-0.5, 0.5)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-2.5, 2.5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-0.5,0.5)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_highSamp_medAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_highSamp_medAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

######################################
######################################
####                              ####
#### High Sample size, High AUC   ####
####                              ####
######################################
######################################

Samp_Cat = "High"
AUC_Cat = "High"
j = 1
# nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",j, ".rds")
dat = readRDS(nms)

True.AUC.BN.equal = dat$fit.BN.equal$TrueAUC
True.AUC.BN.unequal = dat$fit.BN.unequal$TrueAUC
True.AUC.Skewed.I = dat$fit.Skewed.I$TrueAUC
True.AUC.Skewed.II = dat$fit.Skewed.II$TrueAUC
True.AUC.Skewed.III = dat$fit.Skewed.III$TrueAUC
True.AUC.Mixed.I = dat$fit.Mixed.I$TrueAUC
True.AUC.Mixed.II = dat$fit.Mixed.II$TrueAUC

True.ROC.BN.equal = dat$fit.BN.equal$TrueROC
True.ROC.BN.unequal = dat$fit.BN.unequal$TrueROC
True.ROC.Skewed.I = dat$fit.Skewed.I$TrueROC
True.ROC.Skewed.II = dat$fit.Skewed.II$TrueROC
True.ROC.Skewed.III = dat$fit.Skewed.III$TrueROC
True.ROC.Mixed.I = dat$fit.Mixed.I$TrueROC
True.ROC.Mixed.II = dat$fit.Mixed.II$TrueROC

True.Cutoff.BN.equal = dat$fit.BN.equal$True.cutoff
True.Cutoff.BN.unequal = dat$fit.BN.unequal$True.cutoff
True.Cutoff.Skewed.I = dat$fit.Skewed.I$True.cutoff
True.Cutoff.Skewed.II = dat$fit.Skewed.II$True.cutoff
True.Cutoff.Skewed.III = dat$fit.Skewed.III$True.cutoff
True.Cutoff.Mixed.I = dat$fit.Mixed.I$True.cutoff
True.Cutoff.Mixed.II = dat$fit.Mixed.II$True.cutoff

grid = dat$fit.BN.equal$grid
G = 1000

AUC.BN.equal = AUC.BN.unequal = AUC.Skewed.I = AUC.Skewed.II = 
  AUC.Skewed.III = AUC.Mixed.I = AUC.Mixed.II = array(NA, dim = c(G, 7))

ROC.BN.equal = ROC.BN.unequal = ROC.Skewed.I = ROC.Skewed.II = 
  ROC.Skewed.III = ROC.Mixed.I = ROC.Mixed.II = array(NA, dim = c(G, length(grid), 7))

J.BN.equal = J.BN.unequal = J.Skewed.I = J.Skewed.II = 
  J.Skewed.III = J.Mixed.I = J.Mixed.II = array(NA, dim = c(G, 6, 7))
ER.BN.equal = ER.BN.unequal = ER.Skewed.I = ER.Skewed.II = 
  ER.Skewed.III = ER.Mixed.I = ER.Mixed.II = array(NA, dim = c(G, 6, 7))
CZ.BN.equal = CZ.BN.unequal = CZ.Skewed.I = CZ.Skewed.II = 
  CZ.Skewed.III = CZ.Mixed.I = CZ.Mixed.II = array(NA, dim = c(G, 6, 7))
IU.BN.equal = IU.BN.unequal = IU.Skewed.I = IU.Skewed.II = 
  IU.Skewed.III = IU.Mixed.I = IU.Mixed.II = array(NA, dim = c(G, 6, 7))

tt = Sys.time()
for(k in 1:1000){
  
  # nms = paste0("Sim Res/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  nms = paste0("/Users/soutikghosal/Desktop/Simulations/Cutoff/NoCov/res_", Samp_Cat, "Samp_", AUC_Cat, "AUC_",k, ".rds")
  
  
  if(file.exists(nms)){
    dat = readRDS(nms)
    
    AUC.BN.equal[k,] = as.vector(unlist(dat$fit.BN.equal$AUC))
    AUC.BN.unequal[k,] = as.vector(unlist(dat$fit.BN.unequal$AUC))
    AUC.Skewed.I[k,] = as.vector(unlist(dat$fit.Skewed.I$AUC))
    AUC.Skewed.II[k,] = as.vector(unlist(dat$fit.Skewed.II$AUC))
    AUC.Skewed.III[k,] = as.vector(unlist(dat$fit.Skewed.III$AUC))
    AUC.Mixed.I[k,] = as.vector(unlist(dat$fit.Mixed.I$AUC))
    AUC.Mixed.II[k,] = as.vector(unlist(dat$fit.Mixed.II$AUC))
    
    ROC.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ROC[,-1])
    ROC.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ROC[,-1])
    ROC.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ROC[,-1])
    ROC.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ROC[,-1])
    ROC.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ROC[,-1])
    ROC.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ROC[,-1])
    ROC.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ROC[,-1])
    
    J.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$J.tab)
    J.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$J.tab)
    J.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$J.tab)
    J.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$J.tab)
    J.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$J.tab)
    J.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$J.tab)
    J.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$J.tab)
    
    ER.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$ER.tab)
    ER.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$ER.tab)
    ER.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$ER.tab)
    ER.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$ER.tab)
    ER.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$ER.tab)
    ER.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$ER.tab)
    ER.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$ER.tab)
    
    CZ.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$CZ.tab)
    CZ.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$CZ.tab)
    CZ.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$CZ.tab)
    CZ.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$CZ.tab)
    CZ.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$CZ.tab)
    CZ.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$CZ.tab)
    CZ.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$CZ.tab)
    
    IU.BN.equal[k,,] = do.call(cbind, dat$fit.BN.equal$IU.tab)
    IU.BN.unequal[k,,] = do.call(cbind, dat$fit.BN.unequal$IU.tab)
    IU.Skewed.I[k,,] = do.call(cbind, dat$fit.Skewed.I$IU.tab)
    IU.Skewed.II[k,,] = do.call(cbind, dat$fit.Skewed.II$IU.tab)
    IU.Skewed.III[k,,] = do.call(cbind, dat$fit.Skewed.III$IU.tab)
    IU.Mixed.I[k,,] = do.call(cbind, dat$fit.Mixed.I$IU.tab)
    IU.Mixed.II[k,,] = do.call(cbind, dat$fit.Mixed.II$IU.tab)
  } 
}
Sys.time() - tt

mean.NA = function(x){mean(x, na.rm = TRUE)}
sd.NA = function(x){sd(x, na.rm = TRUE)}

BN.equal.summ = summary.tab(AUC = AUC.BN.equal, AUC.True = True.AUC.BN.equal, ROC = ROC.BN.equal,
                            ROC.True = True.ROC.BN.equal, 
                            J = J.BN.equal, ER = ER.BN.equal, 
                            CZ = CZ.BN.equal, IU = IU.BN.equal,
                            Cutoff.tab = True.Cutoff.BN.equal)


BN.unequal.summ = summary.tab(AUC = AUC.BN.unequal, AUC.True = True.AUC.BN.unequal, ROC = ROC.BN.unequal,
                              ROC.True = True.ROC.BN.unequal, 
                              J = J.BN.unequal, ER = ER.BN.unequal, 
                              CZ = CZ.BN.unequal, IU = IU.BN.unequal,
                              Cutoff.tab = True.Cutoff.BN.unequal)

Skewed.I.summ = summary.tab(AUC = AUC.Skewed.I, AUC.True = True.AUC.Skewed.I, ROC = ROC.Skewed.I,
                            ROC.True = True.ROC.Skewed.I, 
                            J = J.Skewed.I, ER = ER.Skewed.I, 
                            CZ = CZ.Skewed.I, IU = IU.Skewed.I,
                            Cutoff.tab = True.Cutoff.Skewed.I)

Skewed.II.summ = summary.tab(AUC = AUC.Skewed.II, AUC.True = True.AUC.Skewed.II, ROC = ROC.Skewed.II,
                             ROC.True = True.ROC.Skewed.II, 
                             J = J.Skewed.II, ER = ER.Skewed.II, 
                             CZ = CZ.Skewed.II, IU = IU.Skewed.II,
                             Cutoff.tab = True.Cutoff.Skewed.II)

Skewed.III.summ = summary.tab(AUC = AUC.Skewed.III, AUC.True = True.AUC.Skewed.III, ROC = ROC.Skewed.III,
                              ROC.True = True.ROC.Skewed.III, 
                              J = J.Skewed.III, ER = ER.Skewed.III, 
                              CZ = CZ.Skewed.III, IU = IU.Skewed.III,
                              Cutoff.tab = True.Cutoff.Skewed.III)


Mixed.I.summ = summary.tab(AUC = AUC.Mixed.I, AUC.True = True.AUC.Mixed.I, ROC = ROC.Mixed.I,
                           ROC.True = True.ROC.Mixed.I, 
                           J = J.Mixed.I, ER = ER.Mixed.I, 
                           CZ = CZ.Mixed.I, IU = IU.Mixed.I,
                           Cutoff.tab = True.Cutoff.Mixed.I)

Mixed.II.summ = summary.tab(AUC = AUC.Mixed.II, AUC.True = True.AUC.Mixed.II, ROC = ROC.Mixed.II,
                            ROC.True = True.ROC.Mixed.II, 
                            J = J.Mixed.II, ER = ER.Mixed.II, 
                            CZ = CZ.Mixed.II, IU = IU.Mixed.II,
                            Cutoff.tab = True.Cutoff.Mixed.II)


bias.tab.me = rbind(bias.table.me(BN.equal.summ),
                    bias.table.me(BN.unequal.summ),
                    bias.table.me(Skewed.I.summ),
                    bias.table.me(Skewed.II.summ),
                    bias.table.me(Skewed.III.summ),
                    bias.table.me(Mixed.I.summ),
                    bias.table.me(Mixed.II.summ)
)

bias.tab.me$Data.gen = c(rep("BN equal", 7),
                         rep("BN unequal", 7),
                         rep("Skewed I", 7),
                         rep("Skewed II", 7),
                         rep("Skewed III", 7),
                         rep("Mixed I", 7),
                         rep("Mixed II", 7))

# write.csv(bias.tab.me, "bias_table_me_HighSamp_HighAUC.csv")

require(ggpubr)
allbutAUC.noCon.bias.plot = ggarrange(BN.equal.summ$allbutAUC.noCon.plot + ylim(c(-0.25, 0.25)),
                                      BN.unequal.summ$allbutAUC.noCon.plot + ylim(c(-0.5, 0.5)),
                                      Skewed.I.summ$allbutAUC.noCon.plot + ylim(c(-7.5, 7.5)),
                                      Skewed.II.summ$allbutAUC.noCon.plot + ylim(c(-5, 5)),
                                      Skewed.III.summ$allbutAUC.noCon.plot + ylim(c(-7.5, 7.5)),
                                      Mixed.I.summ$allbutAUC.noCon.plot + ylim(c(-1,1)),
                                      Mixed.II.summ$allbutAUC.noCon.plot + ylim(c(-0.5,0.5)), nrow = 4, ncol = 2,
                                      common.legend = TRUE,
                                      labels = LETTERS[1:7])

AUC.noCon.bias.plot = ggarrange(BN.equal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                BN.unequal.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Skewed.III.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.I.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                Mixed.II.summ$AUC.noCon.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4),
                                nrow = 4, ncol = 2,
                                common.legend = TRUE,
                                labels = LETTERS[1:7])

ggsave(filename = "Plot_sim_bias_noCon_highSamp_highAUC.png",
       plot = allbutAUC.noCon.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "Plot_sim_AUCbias_noCon_highSamp_highAUC.png",
       plot = AUC.noCon.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)
