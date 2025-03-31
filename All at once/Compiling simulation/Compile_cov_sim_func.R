
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

#################################################
####
#### Choice of sample size (Low, Medium, High)
####
#################################################

Samp_Cat = "Medium"  ## Sample size

comp.cov.sim.func = function(Samp_Cat = "Medium"){
  
  mean.NA = function(x){mean(x, na.rm = TRUE)}
  sd.NA = function(x){sd(x, na.rm = TRUE)}
  
  summary.tab = function(AUC, AUC.True, ROC, ROC.True, 
                         J, ER, CZ, IU, Cutoff.tab){
    
    # AUC = AUC.BN
    # ROC = ROC.BN
    # AUC.True = True.AUC.BN
    # ROC.True = True.ROC.BN
    AUC.hat = apply(AUC,c(2,3),mean.NA)
    AUC.sd = apply(AUC,c(2,3),sd.NA)
    AUC.bias = apply(AUC-AUC.True, c(2,3),mean.NA)
    AUC.bias.SD = apply(AUC-AUC.True,c(2,3),sd.NA)
    AUC.mse = apply((AUC-AUC.True)^2, c(2,3),mean.NA)
    
    EMSE1 = apply(apply((ROC[,,,1]-ROC.True[,1])^2,c(2,3),mean.NA),2,mean.NA)
    EMSE2 = apply(apply((ROC[,,,2]-ROC.True[,2])^2,c(2,3),mean.NA),2,mean.NA)
    
    AUC0.tab = data.frame(True.AUC = AUC.True[1],
                          AUC = AUC.hat[,1], SD = AUC.sd[,1], Bias = AUC.bias[,1], Bias.SD = AUC.bias.SD[,1],
                          MSE = AUC.mse[,1], EMSE = EMSE1)
    AUC1.tab = data.frame(True.AUC = AUC.True[2],
                          AUC = AUC.hat[,2], SD = AUC.sd[,2], Bias = AUC.bias[,2], Bias.SD = AUC.bias.SD[,2],
                          MSE = AUC.mse[,2], EMSE = EMSE2)
    rownames(AUC0.tab)= rownames(AUC1.tab) = c("BN","PV", "Semi.PV")
    AUC.tab = list(AUC0.tab,AUC1.tab)
    
    ROC.hat = apply(ROC,c(2,3,4),mean.NA)
    colnames(ROC.hat) = c("BN","PV", "Semi.PV")
    
    ROC.hat = lapply(seq(dim(ROC.hat)[3]), function(x) ROC.hat[ , , x])
    
    AUC1.bias.whole = AUC[,,1]-AUC.True[1]
    AUC2.bias.whole = AUC[,,2]-AUC.True[2]
    colnames(AUC1.bias.whole) = colnames(AUC2.bias.whole) = c("BN","PV", "Semi.PV")
    AUC.bias.whole = list(AUC1.bias.whole,AUC2.bias.whole)
    
    # J = J.BN
    # ER = ER.BN
    # CZ = CZ.BN
    # IU = IU.BN
    # Cutoff.tab = True.Cutoff.BN
    
    J.hat = list(apply(J[,,,1], c(2,3), mean.NA), apply(J[,,,2], c(2,3), mean.NA))
    ER.hat = list(apply(ER[,,,1], c(2,3), mean.NA), apply(ER[,,,2], c(2,3), mean.NA))
    CZ.hat = list(apply(CZ[,,,1], c(2,3), mean.NA), apply(CZ[,,,2], c(2,3), mean.NA))
    IU.hat = list(apply(IU[,,,1], c(2,3), mean.NA), apply(IU[,,,2], c(2,3), mean.NA))
    
    Bias.J.hat = list(apply(J[,,,1] - Cutoff.tab$cutoff0[1,], c(2,3), mean.NA), 
                      apply(J[,,,2] - Cutoff.tab$cutoff1[1,], c(2,3), mean.NA))
    Bias.ER.hat = list(apply(ER[,,,1] - Cutoff.tab$cutoff0[2,], c(2,3), mean.NA), 
                       apply(ER[,,,2] - Cutoff.tab$cutoff1[2,], c(2,3), mean.NA))
    Bias.CZ.hat = list(apply(CZ[,,,1] - Cutoff.tab$cutoff0[3,], c(2,3), mean.NA), 
                       apply(CZ[,,,2] - Cutoff.tab$cutoff1[3,], c(2,3), mean.NA))
    Bias.IU.hat = list(apply(IU[,,,1] - Cutoff.tab$cutoff0[4,], c(2,3), mean.NA), 
                       apply(IU[,,,2] - Cutoff.tab$cutoff1[4,], c(2,3), mean.NA))
    
    Bias.J.SD = list(apply(J[,,,1] - Cutoff.tab$cutoff0[1,], c(2,3), sd.NA), 
                     apply(J[,,,2] - Cutoff.tab$cutoff1[1,], c(2,3), sd.NA))
    Bias.ER.SD = list(apply(ER[,,,1] - Cutoff.tab$cutoff0[2,], c(2,3), sd.NA), 
                      apply(ER[,,,2] - Cutoff.tab$cutoff1[2,], c(2,3), sd.NA))
    Bias.CZ.SD = list(apply(CZ[,,,1] - Cutoff.tab$cutoff0[3,], c(2,3), sd.NA), 
                      apply(CZ[,,,2] - Cutoff.tab$cutoff1[3,], c(2,3), sd.NA))
    Bias.IU.SD = list(apply(IU[,,,1] - Cutoff.tab$cutoff0[4,], c(2,3), sd.NA), 
                      apply(IU[,,,2] - Cutoff.tab$cutoff1[4,], c(2,3), sd.NA))
    
    
    J.cutoff.bias = list(J[,1,,1] - Cutoff.tab$cutoff0[1,1], J[,1,,2] - Cutoff.tab$cutoff1[1,1])
    ER.cutoff.bias = list(ER[,1,,1] - Cutoff.tab$cutoff0[2,1], ER[,1,,2] - Cutoff.tab$cutoff1[2,1])
    CZ.cutoff.bias = list(CZ[,1,,1] - Cutoff.tab$cutoff0[3,1], CZ[,1,,2] - Cutoff.tab$cutoff1[3,1])
    IU.cutoff.bias = list(IU[,1,,1] - Cutoff.tab$cutoff0[4,1], IU[,1,,2] - Cutoff.tab$cutoff1[4,1])
    
    # colnames(AUC1.bias.whole) = colnames(AUC2.bias.whole)
    colnames(J.cutoff.bias[[1]]) = colnames(J.cutoff.bias[[2]]) = 
      colnames(ER.cutoff.bias[[1]]) = colnames(ER.cutoff.bias[[2]]) = 
      colnames(CZ.cutoff.bias[[1]]) = colnames(CZ.cutoff.bias[[2]]) =
      colnames(IU.cutoff.bias[[1]]) = colnames(IU.cutoff.bias[[2]]) = c("BN","PV", "Semi.PV")
    
    AUC0.bias.2 = melt(AUC1.bias.whole)
    AUC1.bias.2 = melt(AUC2.bias.whole)
    J0.cutoff.bias.2 = melt(J.cutoff.bias[[1]])
    J1.cutoff.bias.2 = melt(J.cutoff.bias[[2]])
    ER0.cutoff.bias.2 = melt(ER.cutoff.bias[[1]])
    ER1.cutoff.bias.2 = melt(ER.cutoff.bias[[2]])
    CZ0.cutoff.bias.2 = melt(CZ.cutoff.bias[[1]])
    CZ1.cutoff.bias.2 = melt(CZ.cutoff.bias[[2]])
    IU0.cutoff.bias.2 = melt(IU.cutoff.bias[[1]])
    IU1.cutoff.bias.2 = melt(IU.cutoff.bias[[2]])
    AUC0.bias.2$Method = AUC1.bias.2$Method = "AUC"
    J0.cutoff.bias.2$Method = J1.cutoff.bias.2$Method = "J"
    ER0.cutoff.bias.2$Method = ER1.cutoff.bias.2$Method = "ER"
    CZ0.cutoff.bias.2$Method = CZ1.cutoff.bias.2$Method = "CZ"
    IU0.cutoff.bias.2$Method = IU1.cutoff.bias.2$Method = "IU"
    
    bias0.whole = rbind(AUC0.bias.2, J0.cutoff.bias.2, ER0.cutoff.bias.2,
                        CZ0.cutoff.bias.2, IU0.cutoff.bias.2)
    bias1.whole = rbind(AUC1.bias.2, J1.cutoff.bias.2, ER1.cutoff.bias.2,
                        CZ1.cutoff.bias.2, IU1.cutoff.bias.2)
    colnames(bias0.whole)[2:3] = colnames(bias1.whole)[2:3] = c("Model", "Bias")
    bias0.whole$X1 = bias1.whole$X1 = NULL
    
    for(i in seq_along(J.hat)) {
      colnames(J.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(ER.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(CZ.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(IU.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(Bias.J.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(Bias.ER.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(Bias.CZ.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      colnames(Bias.IU.hat[[i]]) <- c("BN", "PV", "Semi.PV")
      
      rownames(J.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                "sensitivity (data)", "specificity (data)")
      rownames(ER.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                 "sensitivity (data)", "specificity (data)")
      rownames(CZ.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                 "sensitivity (data)", "specificity (data)")
      rownames(IU.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                 "sensitivity (data)", "specificity (data)")
      rownames(Bias.J.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                     "sensitivity (data)", "specificity (data)")
      rownames(Bias.ER.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                      "sensitivity (data)", "specificity (data)")
      rownames(Bias.CZ.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                      "sensitivity (data)", "specificity (data)")
      rownames(Bias.IU.hat[[i]]) <- c("cutoff","value","sensitivity", "specificity", 
                                      "sensitivity (data)", "specificity (data)")
    }
    
    line_size <- 1 # defining variable upfront as we will re-use it
    base_size <- 16 # defining separately, same as for line_size
    axis_text_rel_size = -1
    title_text_rel_size = +2
    
    # AUC0.plot = ggplot(bias0.whole %>% filter(Method == "AUC"), aes(x = Model, y=Bias)) +
    #   geom_boxplot()
    J0.plot = ggplot(bias0.whole %>% filter(Method == "J"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    ER0.plot = ggplot(bias0.whole %>% filter(Method == "ER"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    CZ0.plot = ggplot(bias0.whole %>% filter(Method == "CZ"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    IU0.plot = ggplot(bias0.whole %>% filter(Method == "IU"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    
    # AUC1.plot = ggplot(bias1.whole %>% filter(Method == "AUC"), aes(x = Model, y=Bias)) +
    #   geom_boxplot()
    J1.plot = ggplot(bias1.whole %>% filter(Method == "J"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    ER1.plot = ggplot(bias1.whole %>% filter(Method == "ER"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    CZ1.plot = ggplot(bias1.whole %>% filter(Method == "CZ"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    IU1.plot = ggplot(bias1.whole %>% filter(Method == "IU"), aes(x = Model, y=Bias)) +
      geom_boxplot()
    
    whole0.plot = ggplot(bias0.whole, aes(x = Model, y=Bias, fill = Method)) +
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
        # axis.title.y = element_text(angle = 90, vjust = 2),
        # axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size=30),
        plot.background = element_blank())
    
    whole1.plot = ggplot(bias1.whole, aes(x = Model, y=Bias, fill = Method)) +
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
        # axis.title.y = element_text(angle = 90, vjust = 2),
        # axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size=30),
        plot.background = element_blank())
    
    allbutAUC0.plot = ggplot(bias0.whole %>% filter(Method != "AUC"), aes(x = Model, y=Bias, fill = Method)) +
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
    
    allbutAUC1.plot = ggplot(bias1.whole %>% filter(Method != "AUC"), aes(x = Model, y=Bias, fill = Method)) +
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
    
    
    AUC0.plot = ggplot(bias0.whole %>% filter(Method == "AUC"), aes(x = Model, y=Bias)) +
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
                                  size = rel((title_text_rel_size + base_size) / (2*base_size)), hjust = 0.5),
        axis.line = element_line(colour="black", linewidth = line_size),
        axis.ticks = element_line(colour="black", linewidth = line_size),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / (2*base_size))),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.background = element_blank())
    
    AUC1.plot = ggplot(bias1.whole %>% filter(Method == "AUC"), aes(x = Model, y=Bias)) +
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
                                  size = rel((title_text_rel_size + base_size) / (2*base_size)), hjust = 0.5),
        axis.line = element_line(colour="black", linewidth = line_size),
        axis.ticks = element_line(colour="black", linewidth = line_size),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / (2*base_size))),
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
               bias0.whole = bias0.whole, bias1.whole = bias1.whole,
               AUC0.plot=AUC0.plot, J0.plot=J0.plot, ER0.plot=ER0.plot, CZ0.plot=CZ0.plot, 
               IU0.plot=IU0.plot,
               AUC1.plot=AUC1.plot, J1.plot=J1.plot, ER1.plot=ER1.plot, CZ1.plot=CZ1.plot, 
               IU1.plot=IU1.plot,
               whole0.plot = whole0.plot, whole1.plot = whole1.plot,
               allbutAUC0.plot = allbutAUC0.plot, allbutAUC1.plot = allbutAUC1.plot)
    return(out)
  }
  
  bias.me.table = function(summ.object){
    
    library(dplyr)
    median.NA = function(x){as.numeric(quantile(x, probs = c(0.5), na.rm = TRUE))}
    IQR.NA = function(x){IQR(x, na.rm = TRUE, type = 7)}
    
    # summ.object = BN.summ
    
    nms = colnames(summ.object$AUC.bias.whole[[1]])
    
    # AUC1.bias = paste(round(AUC1.bias$Bias,3), " ± ", round(AUC1.bias$Bias.IQR,3), sep = "")
    AUC1.bias = paste(round(apply(summ.object$AUC.bias.whole[[1]],2,median.NA),3), " ± ", round(apply(summ.object$AUC.bias.whole[[1]],2,IQR.NA),3), sep = "")
    # AUC2.bias = paste(round(AUC2.bias$Bias,3), " ± ", round(AUC2.bias$Bias.IQR,3), sep = "")
    AUC2.bias = paste(round(apply(summ.object$AUC.bias.whole[[2]],2,median.NA),3), " ± ", round(apply(summ.object$AUC.bias.whole[[2]],2,IQR.NA),3), sep = "")
    
    # J1.Bias = data.frame(Bias = summ.object$Bias.J[[1]][1,], Bias.IQR = summ.object$Bias.J.SD[[1]][1,])
    J1.Bias = data.frame(Bias = apply(summ.object$J.cutoff.bias[[1]],2,median.NA), Bias.IQR = apply(summ.object$J.cutoff.bias[[1]],2,IQR.NA))
    J1.Bias = paste(round(J1.Bias$Bias,3), " ± ", round(J1.Bias$Bias.IQR,3), sep = "")
    # J2.Bias = data.frame(Bias = summ.object$Bias.J[[2]][1,], Bias.IQR = summ.object$Bias.J.SD[[2]][1,])
    J2.Bias = data.frame(Bias = apply(summ.object$J.cutoff.bias[[2]],2,median.NA), Bias.IQR = apply(summ.object$J.cutoff.bias[[2]],2,IQR.NA))
    J2.Bias = paste(round(J2.Bias$Bias,3), " ± ", round(J2.Bias$Bias.IQR,3), sep = "")
    
    ER1.Bias = data.frame(Bias = apply(summ.object$ER.cutoff.bias[[1]],2,median.NA), Bias.IQR = apply(summ.object$ER.cutoff.bias[[1]],2,IQR.NA))
    ER1.Bias = paste(round(ER1.Bias$Bias,3), " ± ", round(ER1.Bias$Bias.IQR,3), sep = "")
    ER2.Bias = data.frame(Bias = apply(summ.object$ER.cutoff.bias[[2]],2,median.NA), Bias.IQR = apply(summ.object$ER.cutoff.bias[[2]],2,IQR.NA))
    ER2.Bias = paste(round(ER2.Bias$Bias,3), " ± ", round(ER2.Bias$Bias.IQR,3), sep = "")
    
    CZ1.Bias = data.frame(Bias = apply(summ.object$CZ.cutoff.bias[[1]],2,median.NA), Bias.IQR = apply(summ.object$CZ.cutoff.bias[[1]],2,IQR.NA))
    CZ1.Bias = paste(round(CZ1.Bias$Bias,3), " ± ", round(CZ1.Bias$Bias.IQR,3), sep = "")
    CZ2.Bias = data.frame(Bias = apply(summ.object$CZ.cutoff.bias[[2]],2,median.NA), Bias.IQR = apply(summ.object$CZ.cutoff.bias[[2]],2,IQR.NA))
    CZ2.Bias = paste(round(CZ2.Bias$Bias,3), " ± ", round(CZ2.Bias$Bias.IQR,3), sep = "")
    
    IU1.Bias = data.frame(Bias = apply(summ.object$IU.cutoff.bias[[1]],2,median.NA), Bias.IQR = apply(summ.object$IU.cutoff.bias[[1]],2,IQR.NA))
    IU1.Bias = paste(round(IU1.Bias$Bias,3), " ± ", round(IU1.Bias$Bias.IQR,3), sep = "")
    IU2.Bias = data.frame(Bias = apply(summ.object$IU.cutoff.bias[[2]],2,median.NA), Bias.IQR = apply(summ.object$IU.cutoff.bias[[2]],2,IQR.NA))
    IU2.Bias = paste(round(IU2.Bias$Bias,3), " ± ", round(IU2.Bias$Bias.IQR,3), sep = "")
    
    final1.bias = data.frame(AUC = AUC1.bias, J = J1.Bias, ER = ER1.Bias,
                             CZ = CZ1.Bias, IU = IU1.Bias)
    final2.bias = data.frame(AUC = AUC2.bias, J = J2.Bias, ER = ER2.Bias,
                             CZ = CZ2.Bias, IU = IU2.Bias)
    rownames(final1.bias) = rownames(final2.bias) = nms
    final.bias = list(final1.bias, final2.bias)
    return(final.bias)
  }
  
  fit.name = paste0("SimData_cov_",Samp_Cat, "Samp.rds")
  fit = readRDS(fit.name)
  
  BN.summ = summary.tab(AUC = fit$AUC.est$AUC.BN, 
                        AUC.True = fit$True.AUC$True.AUC.BN, 
                        ROC = fit$ROC.est$ROC.BN,
                        ROC.True = fit$True.ROC$True.ROC.BN, 
                        J = fit$J.est$J.BN, ER =fit$ER.est$ER.BN, 
                        CZ = fit$CZ.est$CZ.BN, fit$IU.est$IU.BN,
                        Cutoff.tab = fit$True.Cutoff$True.Cutoff.BN)
  
  Skewed.summ = summary.tab(AUC = fit$AUC.est$AUC.Skewed, 
                            AUC.True = fit$True.AUC$True.AUC.Skewed, 
                            ROC = fit$ROC.est$ROC.Skewed,
                            ROC.True = fit$True.ROC$True.ROC.Skewed, 
                            J = fit$J.est$J.Skewed, ER =fit$ER.est$ER.Skewed, 
                            CZ = fit$CZ.est$CZ.Skewed, fit$IU.est$IU.Skewed,
                            Cutoff.tab = fit$True.Cutoff$True.Cutoff.Skewed)
  
  Mixed.summ = summary.tab(AUC = fit$AUC.est$AUC.Mixed, 
                           AUC.True = fit$True.AUC$True.AUC.Mixed, 
                           ROC = fit$ROC.est$ROC.Mixed,
                           ROC.True = fit$True.ROC$True.ROC.Mixed, 
                           J = fit$J.est$J.Mixed, ER =fit$ER.est$ER.Mixed, 
                           CZ = fit$CZ.est$CZ.Mixed, fit$IU.est$IU.Mixed,
                           Cutoff.tab = fit$True.Cutoff$True.Cutoff.Mixed)
  
  
  bias.tab.me.x0 = rbind(bias.me.table(BN.summ)[[1]],
                         bias.me.table(Skewed.summ)[[1]],
                         bias.me.table(Mixed.summ)[[1]]
  )
  
  bias.tab.me.x1 = rbind(bias.me.table(BN.summ)[[2]],
                         bias.me.table(Skewed.summ)[[2]],
                         bias.me.table(Mixed.summ)[[2]]
  )
  
  bias.me.tab = rbind(bias.tab.me.x0,bias.tab.me.x1)
  bias.me.tab$Data.gen = rep(c(rep("BN", 3),
                           rep("Skewed", 3),
                           rep("Mixed", 3)),2)
  bias.me.tab$Fit.mod = rep(c("BN","PV","Semi.PV"),6)
    
  # tab.name = paste0("bias_table_me_cov_",Samp_Cat,"Samp.csv")
  
  require(ggpubr)
  allbutAUC.bias.plot = ggarrange(BN.summ$allbutAUC0.plot + ylim(c(-1, 1)),
                                  BN.summ$allbutAUC1.plot + ylim(c(-1, 1)),
                                  Skewed.summ$allbutAUC0.plot + ylim(c(-10, 10)),
                                  Skewed.summ$allbutAUC1.plot + ylim(c(-10, 10)),
                                  Mixed.summ$allbutAUC0.plot + ylim(c(-1, 1)),
                                  Mixed.summ$allbutAUC1.plot + ylim(c(-1, 1)), nrow = 3, ncol = 2,
                                  common.legend = TRUE,
                                  labels = LETTERS[1:6])
  
  AUC.bias.plot = ggarrange(BN.summ$AUC0.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) + ylim(c(-0.2, 0.2)),
                            BN.summ$AUC1.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) + ylim(c(-0.2, 0.2)),
                            Skewed.summ$AUC0.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) + ylim(c(-0.2, 0.2)),
                            Skewed.summ$AUC1.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) + ylim(c(-0.2, 0.2)),
                            Mixed.summ$AUC0.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) + ylim(c(-0.2, 0.2)),
                            Mixed.summ$AUC1.plot + geom_hline(yintercept=0, linetype="dashed", alpha = 0.4) + ylim(c(-0.2, 0.2)), 
                            nrow = 3, ncol = 2,
                            common.legend = TRUE,
                            labels = LETTERS[1:6])
  
  # Bias.name = paste0("Plot_sim_bias_cov_",Samp_Cat,"Samp.png")
  # AUC.Bias.name = paste0("Plot_sim_AUCbias_cov_",Samp_Cat,"Samp.png")
  
  out = list(bias.tab.me = bias.me.tab,
             cutoff.bias.plot = allbutAUC.bias.plot,
             AUC.bias.plot = AUC.bias.plot)
  return(out)
}


## The following creates Figure B.19 in the supplement
## Figures B.1 and B.5 can be plotted by assuming
## AUC_Cat = "Low" and "High" respectively
## Other similar plots with different sample sizes
## can be made similarly.

# ggsave(filename = Bias.name,
#        plot = allbutAUC.bias.plot,
#        width = 12, height = 16,
#        device='png' #, dpi=1000
# )

# ggsave(filename = "FigureB19.csv",
#        plot = allbutAUC.bias.plot,
#        width = 12, height = 16,
#        device='png' #, dpi=1000
# )

# 
# ggsave(filename = AUC.Bias.name,
#        plot = AUC.bias.plot,
#        width = 7, height = 8,
#        device='png' #, dpi=1000
# )
# ggsave(filename = "FigureB20.csv",
#        plot = AUC.bias.plot,
#        width = 7, height = 8,
#        device='png' #, dpi=1000
# )
