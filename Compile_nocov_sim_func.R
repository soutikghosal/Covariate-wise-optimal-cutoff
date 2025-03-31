
##########################
######
###### This script contains
###### codes that summarizes
###### the existing simulations 
###### for no covariate framework
######
##########################


rm(list = ls())

library(reshape)
library(ggplot2)
library(dplyr)
library(ggthemes)

#################################################
####
#### Choice of sample size (Low, Medium, High)
#### Choice of AUC level (Low, Medium, High)
####
#################################################

Samp_Cat = "Medium" ## Sample size
AUC_Cat = "Medium" ## AUC level

comp.nocov.sim.func = function(Samp_Cat = "Medium", AUC_Cat = "Medium"){
  
  mean.NA = function(x){mean(x, na.rm = TRUE)}
  sd.NA = function(x){sd(x, na.rm = TRUE)}
  
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
    base_size <- 14 # defining separately, same as for line_size
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
  
  bias.table.me = function(summ.object){
    
    median.NA = function(x){as.numeric(quantile(x, probs = c(0.5), na.rm = TRUE))}
    IQR.NA = function(x){IQR(x, na.rm = TRUE, type = 7)}
    
    library(dplyr)
    # summ.object = BN.unequal.summ
    
    AUC.bias = list(Bias = apply(summ.object$AUC.bias.whole,2,median.NA),
                    Bias.IQR = apply(summ.object$AUC.bias.whole,2,IQR.NA))
    
    # AUC.bias = summ.object$AUC.tab %>% select(Bias, Bias.IQR)
    nms = names(AUC.bias$Bias)
    AUC.bias = paste(round(AUC.bias$Bias,3), " ± ", round(AUC.bias$Bias.IQR,3), sep = "")
    
    # J.Bias = data.frame(Bias = summ.object$Bias.J[1,], Bias.IQR = summ.object$Bias.J.IQR[1,])
    J.Bias = data.frame(Bias = apply(summ.object$J.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$J.cutoff.bias,2,IQR.NA))
    J.Bias = paste(round(J.Bias$Bias,3), " ± ", round(J.Bias$Bias.IQR,3), sep = "")
    
    # ER.Bias = data.frame(Bias = summ.object$Bias.ER[1,], Bias.IQR = summ.object$Bias.ER.IQR[1,])
    ER.Bias = data.frame(Bias = apply(summ.object$ER.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$ER.cutoff.bias,2,IQR.NA))
    ER.Bias = paste(round(ER.Bias$Bias,3), " ± ", round(ER.Bias$Bias.IQR,3), sep = "")
    
    # CZ.Bias = data.frame(Bias = summ.object$Bias.CZ[1,], Bias.IQR = summ.object$Bias.CZ.IQR[1,])
    CZ.Bias = data.frame(Bias = apply(summ.object$CZ.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$CZ.cutoff.bias,2,IQR.NA))
    CZ.Bias = paste(round(CZ.Bias$Bias,3), " ± ", round(CZ.Bias$Bias.IQR,3), sep = "")
    
    # IU.Bias = data.frame(Bias = summ.object$Bias.IU[1,], Bias.IQR = summ.object$Bias.IU.IQR[1,])
    IU.Bias = data.frame(Bias = apply(summ.object$IU.cutoff.bias,2,median.NA), Bias.IQR = apply(summ.object$IU.cutoff.bias,2,IQR.NA))
    IU.Bias = paste(round(IU.Bias$Bias,3), " ± ", round(IU.Bias$Bias.IQR,3), sep = "")
    
    final.bias = data.frame(AUC = AUC.bias, J = J.Bias, ER = ER.Bias,
                            CZ = CZ.Bias, IU = IU.Bias)
    rownames(final.bias) = nms
    return(final.bias)
  }
  
  fit.name = paste0("SimData_",Samp_Cat, "Samp_", AUC_Cat, "AUC.rds")
  fit = readRDS(fit.name)
  
  BN.equal.summ = summary.tab(AUC = fit$AUC.est$AUC.BN.equal, 
                              AUC.True = fit$True.AUC$True.AUC.BN.equal, 
                              ROC = fit$ROC.est$ROC.BN.equal,
                              ROC.True = fit$True.ROC$True.ROC.BN.equal, 
                              J = fit$J.est$J.BN.equal, ER = fit$ER.est$ER.BN.equal, 
                              CZ = fit$CZ.est$CZ.BN.equal, IU = fit$IU.est$IU.BN.equal,
                              Cutoff.tab = fit$True.Cutoff$True.Cutoff.BN.equal)
  
  BN.unequal.summ = summary.tab(AUC = fit$AUC.est$AUC.BN.unequal, 
                                AUC.True = fit$True.AUC$True.AUC.BN.unequal, 
                                ROC = fit$ROC.est$ROC.BN.unequal,
                                ROC.True = fit$True.ROC$True.ROC.BN.unequal, 
                                J = fit$J.est$J.BN.unequal, ER = fit$ER.est$ER.BN.unequal, 
                                CZ = fit$CZ.est$CZ.BN.unequal, IU = fit$IU.est$IU.BN.unequal,
                                Cutoff.tab = fit$True.Cutoff$True.Cutoff.BN.unequal)
  
  Skewed.I.summ = summary.tab(AUC = fit$AUC.est$AUC.Skewed.I, 
                              AUC.True = fit$True.AUC$True.AUC.Skewed.I, 
                              ROC = fit$ROC.est$ROC.Skewed.I,
                              ROC.True = fit$True.ROC$True.ROC.Skewed.I, 
                              J = fit$J.est$J.Skewed.I, ER = fit$ER.est$ER.Skewed.I, 
                              CZ = fit$CZ.est$CZ.Skewed.I, IU = fit$IU.est$IU.Skewed.I,
                              Cutoff.tab = fit$True.Cutoff$True.Cutoff.Skewed.I)
  
  Skewed.II.summ = summary.tab(AUC = fit$AUC.est$AUC.Skewed.II, 
                               AUC.True = fit$True.AUC$True.AUC.Skewed.II, 
                               ROC = fit$ROC.est$ROC.Skewed.II,
                               ROC.True = fit$True.ROC$True.ROC.Skewed.II, 
                               J = fit$J.est$J.Skewed.II, ER = fit$ER.est$ER.Skewed.II, 
                               CZ = fit$CZ.est$CZ.Skewed.II, IU = fit$IU.est$IU.Skewed.II,
                               Cutoff.tab = fit$True.Cutoff$True.Cutoff.Skewed.II)
  
  Skewed.III.summ = summary.tab(AUC = fit$AUC.est$AUC.Skewed.III, 
                                AUC.True = fit$True.AUC$True.AUC.Skewed.III, 
                                ROC = fit$ROC.est$ROC.Skewed.III,
                                ROC.True = fit$True.ROC$True.ROC.Skewed.III, 
                                J = fit$J.est$J.Skewed.III, ER = fit$ER.est$ER.Skewed.III, 
                                CZ = fit$CZ.est$CZ.Skewed.III, IU = fit$IU.est$IU.Skewed.III,
                                Cutoff.tab = fit$True.Cutoff$True.Cutoff.Skewed.III)
  
  Mixed.I.summ = summary.tab(AUC = fit$AUC.est$AUC.Mixed.I, 
                             AUC.True = fit$True.AUC$True.AUC.Mixed.I, 
                             ROC = fit$ROC.est$ROC.Mixed.I,
                             ROC.True = fit$True.ROC$True.ROC.Mixed.I, 
                             J = fit$J.est$J.Mixed.I, ER = fit$ER.est$ER.Mixed.I, 
                             CZ = fit$CZ.est$CZ.Mixed.I, IU = fit$IU.est$IU.Mixed.I,
                             Cutoff.tab = fit$True.Cutoff$True.Cutoff.Mixed.I)
  
  Mixed.II.summ = summary.tab(AUC = fit$AUC.est$AUC.Mixed.II, 
                              AUC.True = fit$True.AUC$True.AUC.Mixed.II, 
                              ROC = fit$ROC.est$ROC.Mixed.II,
                              ROC.True = fit$True.ROC$True.ROC.Mixed.II, 
                              J = fit$J.est$J.Mixed.II, ER = fit$ER.est$ER.Mixed.II, 
                              CZ = fit$CZ.est$CZ.Mixed.II, IU = fit$IU.est$IU.Mixed.II,
                              Cutoff.tab = fit$True.Cutoff$True.Cutoff.Mixed.II)
  
  
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
  
  ## rows to remove as it is
  ## not part of the simulation
  
  BG.id = seq(3,45, length = 7) ## Previously fit BG model
  BiChi.id = BG.id +1 ## Previously fit Proper BN model
  all.id = sort(c(BG.id, BiChi.id)) ## remove these results
  bias.tab.me = bias.tab.me[-all.id,]
  
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
  
  # Bias.name = paste0("Plot_sim_bias_",Samp_Cat,"Samp_", AUC_Cat,"AUC.png")
  # AUC.Bias.name = paste0("Plot_sim_AUCbias_",Samp_Cat,"Samp_", AUC_Cat,"AUC.png")
  
  out = list(bias.tab.me = bias.tab.me,
             cutoff.bias.plot = allbutAUC.noCon.bias.plot,
             AUC.bias.plot = AUC.noCon.bias.plot)
  return(out)
}



## This will create Table 3.b
## Similarly, Tables 3.a and 3.c can be obtained 
## by using Samp_Cat = "Medium" and
## AUC_Cat = "Low" and "High" respectively

## To obtain Tables B.3a - B.3c in the supplement, 
## we will similarly use Samp_Cat = "Low", and
## To obtain Tables B.4a - B.4c, we will similarly use
## Samp_Cat = "High"


## The following creates Figure B.3 in the supplement
## Figures B.1 and B.5 can be plotted by assuming
## AUC_Cat = "Low" and "High" respectively
## Other similar plots with different sample sizes
## can be made similarly.



## The following creates Figure B.4 in the supplement
## Figures B.2 and B.6 can be plotted by assuming
## AUC_Cat = "Low" and "High" respectively
## Other similar plots with different sample sizes
## can be made similarly.


# ggsave(filename = AUC.Bias.name,
#        plot = AUC.noCon.bias.plot,
#        width = 7, height = 8,
#        device='png' #, dpi=1000
# )

# ggsave(filename = "FigureB4.png",
#        plot = AUC.noCon.bias.plot,
#        width = 7, height = 8,
#        device='png' #, dpi=1000
# )
