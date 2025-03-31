##################################
###General info of the article####
##################################

# Supplementary information / reproducible research files for the manuscript 
# Title: "Impact of methodological assumptions and covariates on the cutoff estimation in ROC analysis"

# Authors: Ghosal, S.
# Code was written by Ghosal, S.
# In case of questions or comments please contact soutik.ghosal@virginia.edu!

# The code was written/evaluated in R with the following software versions:
# R version 4.3.2 (2023-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 15.1.1


#############################################################################################
#############################################################################################
###################                                                       ###################
###################          Part 1: Run Synthetic data analysis          ###################
###################                                                       ###################
#############################################################################################
#############################################################################################

rm(list = ls())
setwd("Synthetic data analysis/")
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
library(abind)

###########################
############################
# No covariate analysis
############################
############################

source("data_func.R")

tt = Sys.time()
source("nocov_data_analysis.R")
Sys.time() - tt


AUC_cutoff_tab
AUC_cutoff_tab = cbind(Biomarker = c(rep("AB42",5),rep("Tau", 5),rep("pTau", 5)), AUC_cutoff_tab)
write.csv(AUC_cutoff_tab, "Table8.csv")

all.ROC.plot

ggsave(filename = "Figure1.png",
       plot = all.ROC.plot,
       width = 12, height = 6,
       device='png' #, dpi=1000
)

all.cutoff.plot

ggsave(filename = "Figure2.png",
       plot = all.cutoff.plot,
       width = 18, height = 12,
       device='png' #, dpi=1000
)


############################
############################
# With covariate analysis
############################
############################

rm(list = ls())

source("data_func_cov.R")

tt = Sys.time()
source("cov_data_analysis.R")
Sys.time() - tt

AUC_cutoff_tab

Covariate.level = rep(c(rep("Sex: Male",3),rep("Sex: Female",3)), 3)
Biomarker = c(rep("AB42",6),rep("Tau",6),rep("pTau",6))
AUC_cutoff_tab = cbind(Biomarker,Covariate.level,AUC_cutoff_tab)
write.csv(AUC_cutoff_tab, "Table9.csv")

all.ROC.plot
ggsave(filename = "Figure3.png",
       plot = all.ROC.plot,
       width = 12, height = 18,
       device='png' #, dpi=1000
)

AB42.whole.summ$all.cutoff.plot
ggsave(filename = "Figure4.png",
       plot = AB42.whole.summ$all.cutoff.plot,
       width = 12, height = 18,
       device='png' #, dpi=1000
)

Tau.whole.summ$all.cutoff.plot
ggsave(filename = "Figure5.png",
       plot = Tau.whole.summ$all.cutoff.plot,
       width = 12, height = 18,
       device='png' #, dpi=1000
)

pTau.whole.summ$all.cutoff.plot
ggsave(filename = "Figure6.png",
       plot = pTau.whole.summ$all.cutoff.plot,
       width = 12, height = 18,
       device='png' #, dpi=1000
)

setwd("../")

#############################################################################################
#############################################################################################
###################                                                       ###################
###################          Part 2: Run Compiled Simulation res          ###################
###################                                                       ###################
#############################################################################################
#############################################################################################

rm(list = ls())
setwd("Compiling simulation/")

############################
############################
# No covariate analysis
############################
############################

rm(list = ls())
library(reshape)
library(ggplot2)
library(dplyr)
library(ggthemes)

source("Compile_nocov_sim_func.R")

#################################################
####
#### Choice of sample size (Low, Medium, High)
#### Choice of AUC level (Low, Medium, High)
####
#################################################


fit.sampLow.AUCLow = comp.nocov.sim.func(Samp_Cat = "Low", AUC_Cat = "Low")
write.csv(fit.sampLow.AUCLow$bias.tab.me, "TableB3a.csv")
ggsave(filename = "FigureB7.png",
       plot = fit.sampLow.AUCLow$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB8.png",
       plot = fit.sampLow.AUCLow$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampLow.AUCMed = comp.nocov.sim.func(Samp_Cat = "Low", AUC_Cat = "Medium")
write.csv(fit.sampLow.AUCMed$bias.tab.me, "TableB3b.csv")
ggsave(filename = "FigureB9.png",
       plot = fit.sampLow.AUCMed$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB10.png",
       plot = fit.sampLow.AUCMed$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampLow.AUCHigh = comp.nocov.sim.func(Samp_Cat = "Low", AUC_Cat = "High")
write.csv(fit.sampLow.AUCHigh$bias.tab.me, "TableB3c.csv")
ggsave(filename = "FigureB11.png",
       plot = fit.sampLow.AUCHigh$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB12.png",
       plot = fit.sampLow.AUCHigh$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)



fit.sampMed.AUCLow = comp.nocov.sim.func(Samp_Cat = "Medium", AUC_Cat = "Low")
write.csv(fit.sampMed.AUCLow$bias.tab.me, "Table3a.csv")
ggsave(filename = "FigureB1.png",
       plot = fit.sampMed.AUCLow$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB2.png",
       plot = fit.sampMed.AUCLow$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampMed.AUCMed = comp.nocov.sim.func(Samp_Cat = "Medium", AUC_Cat = "Medium")
write.csv(fit.sampMed.AUCMed$bias.tab.me, "Table3b.csv")
ggsave(filename = "FigureB3.png",
       plot = fit.sampMed.AUCMed$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB4.png",
       plot = fit.sampMed.AUCMed$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampMed.AUCHigh = comp.nocov.sim.func(Samp_Cat = "Medium", AUC_Cat = "High")
write.csv(fit.sampMed.AUCHigh$bias.tab.me, "Table3c.csv")
ggsave(filename = "FigureB5.png",
       plot = fit.sampMed.AUCHigh$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB6.png",
       plot = fit.sampMed.AUCHigh$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)


fit.sampHigh.AUCLow = comp.nocov.sim.func(Samp_Cat = "High", AUC_Cat = "Low")
write.csv(fit.sampHigh.AUCLow$bias.tab.me, "TableB4a.csv")
ggsave(filename = "FigureB13.png",
       plot = fit.sampHigh.AUCLow$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB14.png",
       plot = fit.sampHigh.AUCLow$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampHigh.AUCMed = comp.nocov.sim.func(Samp_Cat = "High", AUC_Cat = "Medium")
write.csv(fit.sampHigh.AUCMed$bias.tab.me, "TableB4b.csv")
ggsave(filename = "FigureB15.png",
       plot = fit.sampHigh.AUCMed$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB16.png",
       plot = fit.sampHigh.AUCMed$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampHigh.AUCHigh = comp.nocov.sim.func(Samp_Cat = "High", AUC_Cat = "High")
write.csv(fit.sampHigh.AUCHigh$bias.tab.me, "TableB4c.csv")
ggsave(filename = "FigureB17.png",
       plot = fit.sampHigh.AUCHigh$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)
ggsave(filename = "FigureB18.png",
       plot = fit.sampHigh.AUCHigh$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

setwd("../")

############################
############################
# With covariate analysis
############################
############################

rm(list = ls())
setwd("Compiling simulation/")
source("Compile_cov_sim_func.R")
library(reshape)
library(ggplot2)
library(dplyr)

#################################################
####
#### Choice of sample size (Low, Medium, High)
####
#################################################

fit.sampLow = comp.cov.sim.func(Samp_Cat = "Low")

write.csv(fit.sampLow$bias.tab.me, "TableB5.csv")

ggsave(filename = "FigureB21.png",
       plot = fit.sampLow$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "FigureB22.png",
       plot = fit.sampLow$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampMed = comp.cov.sim.func(Samp_Cat = "Medium")

write.csv(fit.sampMed$bias.tab.me, "Table6.csv")

ggsave(filename = "FigureB19.png",
       plot = fit.sampMed$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "FigureB20.png",
       plot = fit.sampMed$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

fit.sampHigh = comp.cov.sim.func(Samp_Cat = "High")

write.csv(fit.sampHigh$bias.tab.me, "TableB6.csv")

ggsave(filename = "FigureB23.png",
       plot = fit.sampHigh$cutoff.bias.plot,
       width = 12, height = 16,
       device='png' #, dpi=1000
)

ggsave(filename = "FigureB24.png",
       plot = fit.sampHigh$AUC.bias.plot,
       width = 7, height = 8,
       device='png' #, dpi=1000
)

setwd("../")

#####################################################################################
#####################################################################################
###################                                               ###################
###################          Part 3: Run New Simulations          ###################
###################                                               ###################
#####################################################################################
#####################################################################################

setwd("Simulation code/")

############################
############################
# No covariate analysis
############################
############################

rm(list = ls())
source("Nocov sim func.R")

## One particular example:
## n = 50
## sep = "Low"

n.iterations = 3 ## Number of simulated data
seeds = 1:n.iterations
n = 50
sep_cat = "Low"

True.AUC = array(NA, dim = c(1,7, n.iterations))
True.ROC = array(NA, dim = c(100, 8, n.iterations))
True.Cutoff = array(NA, dim = c(7,4, n.iterations))

AUC.tab = J.tab = ER.tab = CZ.tab = IU.tab = array(NA, dim = c(7,5, n.iterations))
ROC.tab = array(NA, dim = c(100, 5, 7, n.iterations))
# ROC.tab = c()

tt = Sys.time()
for(k in 1:n.iterations){
  
  ff = final.fit(seed = seeds[k], n = n, sep = sep_cat)
  
  ## True values (not changing at different loop)

  TrueAUC = cbind(ff$fit.BN.equal$TrueAUC,
                  ff$fit.BN.unequal$TrueAUC,
                  ff$fit.Skewed.I$TrueAUC,
                  ff$fit.Skewed.II$TrueAUC,
                  ff$fit.Skewed.III$TrueAUC,
                  ff$fit.Mixed.I$TrueAUC,
                  ff$fit.Mixed.II$TrueAUC)
  colnames(TrueAUC) = c("BN.equal","BN.unequal","Skewed.I","Skewed.II",
                        "Skewed.III", "Mixed.I", "Mixed.II")
  TrueROC = cbind(ff$fit.BN.equal$TrueROC,
                  ff$fit.BN.unequal$TrueROC,
                  ff$fit.Skewed.I$TrueROC,
                  ff$fit.Skewed.II$TrueROC,
                  ff$fit.Skewed.III$TrueROC,
                  ff$fit.Mixed.I$TrueROC,
                  ff$fit.Mixed.II$TrueROC)
  colnames(TrueROC) = c("BN.equal","BN.unequal","Skewed.I","Skewed.II",
                        "Skewed.III", "Mixed.I", "Mixed.II")
  TrueROC = data.frame(grid = seq(0,1,length = 100), TrueROC)

  TrueCutoff = rbind(ff$fit.BN.equal$True.cutoff[,1],
                     ff$fit.BN.unequal$True.cutoff[,1],
                     ff$fit.Skewed.I$True.cutoff[,1],
                     ff$fit.Skewed.II$True.cutoff[,1],
                     ff$fit.Skewed.III$True.cutoff[,1],
                     ff$fit.Mixed.I$True.cutoff[,1],
                     ff$fit.Mixed.II$True.cutoff[,1])
  rownames(TrueCutoff) = c("BN.equal","BN.unequal","Skewed.I","Skewed.II",
                         "Skewed.III", "Mixed.I", "Mixed.II")

  True.AUC[,,k] = TrueAUC
  True.ROC[,,k] = as.matrix(TrueROC)
  True.Cutoff[,,k] = as.matrix(TrueCutoff)
  
  ## Estimates
  
  AUC = rbind(ff$fit.BN.equal$AUC,
              ff$fit.BN.unequal$AUC,
              ff$fit.Skewed.I$AUC,
              ff$fit.Skewed.II$AUC,
              ff$fit.Skewed.III$AUC,
              ff$fit.Mixed.I$AUC,
              ff$fit.Mixed.II$AUC)

  ROC = abind(list(ff$fit.BN.equal$ROC[,-1],
                   ff$fit.BN.unequal$ROC[,-1],
                   ff$fit.Skewed.I$ROC[,-1],
                   ff$fit.Skewed.II$ROC[,-1],
                   ff$fit.Skewed.III$ROC[,-1],
                   ff$fit.Mixed.I$ROC[,-1],
                   ff$fit.Mixed.II$ROC[,-1]), along = 3)
  
  J = rbind(ff$fit.BN.equal$J.tab[1,],
            ff$fit.BN.unequal$J.tab[1,],
            ff$fit.Skewed.I$J.tab[1,],
            ff$fit.Skewed.II$J.tab[1,],
            ff$fit.Skewed.III$J.tab[1,],
            ff$fit.Mixed.I$J.tab[1,],
            ff$fit.Mixed.II$J.tab[1,])
  ER = rbind(ff$fit.BN.equal$ER.tab[1,],
            ff$fit.BN.unequal$ER.tab[1,],
            ff$fit.Skewed.I$ER.tab[1,],
            ff$fit.Skewed.II$ER.tab[1,],
            ff$fit.Skewed.III$ER.tab[1,],
            ff$fit.Mixed.I$ER.tab[1,],
            ff$fit.Mixed.II$ER.tab[1,])
  CZ = rbind(ff$fit.BN.equal$CZ.tab[1,],
             ff$fit.BN.unequal$CZ.tab[1,],
             ff$fit.Skewed.I$CZ.tab[1,],
             ff$fit.Skewed.II$CZ.tab[1,],
             ff$fit.Skewed.III$CZ.tab[1,],
             ff$fit.Mixed.I$CZ.tab[1,],
             ff$fit.Mixed.II$CZ.tab[1,])
  IU = rbind(ff$fit.BN.equal$IU.tab[1,],
             ff$fit.BN.unequal$IU.tab[1,],
             ff$fit.Skewed.I$IU.tab[1,],
             ff$fit.Skewed.II$IU.tab[1,],
             ff$fit.Skewed.III$IU.tab[1,],
             ff$fit.Mixed.I$IU.tab[1,],
             ff$fit.Mixed.II$IU.tab[1,])
  
  rownames(AUC) = rownames(J) = rownames(ER) = rownames(CZ) = rownames(IU) = 
    c("BN.equal","BN.unequal","Skewed.I","Skewed.II",
                  "Skewed.III", "Mixed.I", "Mixed.II")
  names(ROC) = c("BN.equal","BN.unequal","Skewed.I","Skewed.II",
                    "Skewed.III", "Mixed.I", "Mixed.II")
  
  ROC.tab[,,,k] = ROC
  AUC.tab[,,k] = as.matrix(AUC)
  J.tab[,,k] = as.matrix(J)
  ER.tab[,,k] = as.matrix(ER)
  CZ.tab[,,k] = as.matrix(CZ)
  IU.tab[,,k] = as.matrix(IU)
  # ROC.tab = c(ROC.tab, ROC)
  
}
Sys.time() - tt

## True values

TrueAUC = True.AUC[,,1]
TrueROC = True.ROC[,,1]
TrueCutoff = True.Cutoff[,,1]

names(TrueAUC) = colnames(TrueROC)[2:8] = rownames(TrueCutoff) = 
  c("BN.equal","BN.unequal","Skewed.I","Skewed.II",
                      "Skewed.III", "Mixed.I", "Mixed.II")
colnames(TrueCutoff) = c("J","ER","CZ","IU")

## Estimated values

AUC.hat = apply(AUC.tab, c(1,2), mean)
ROC.hat = apply(ROC.tab, c(1,2,3), mean)
J.hat = apply(J.tab, c(1,2), mean)
ER.hat = apply(ER.tab, c(1,2), mean)
CZ.hat = apply(CZ.tab, c(1,2), mean)
IU.hat = apply(IU.tab, c(1,2), mean)


############################
############################
# With covariate analysis
############################
############################

rm(list = ls())
source("Cov sim func.R")

## One particular example:
## n = 50
n.iterations = 3 ## Number of simulated data
seeds = 1:n.iterations
n = 50

AUC0.tab = J0.tab = ER0.tab = CZ0.tab = IU0.tab = array(NA, dim = c(3, 3, n.iterations))
AUC1.tab = J1.tab = ER1.tab = CZ1.tab = IU1.tab = array(NA, dim = c(3, 3, n.iterations))
ROC0.tab = ROC1.tab = array(NA, dim = c(100,3,3,n.iterations))

tt = Sys.time()
for(k in 1:n.iterations){
  
  ff = final.fit(seed = seeds[k], n = n)
  
  ## Estimates
  
  AUC0 = as.matrix(rbind(ff$fit.BN$AUC[1,],
                         ff$fit.Skewed$AUC[1,],
                         ff$fit.Mixed$AUC[1,]))
  AUC1 = as.matrix(rbind(ff$fit.BN$AUC[2,],
                         ff$fit.Skewed$AUC[2,],
                         ff$fit.Mixed$AUC[2,]))
  rownames(AUC0) = rownames(AUC1) = c("BN","Skewed","Mixed")
  
  ROC0 = abind(list(BN = data.frame(BN = ff$fit.BN$ROC$BN[,2],
                              PV = ff$fit.BN$ROC$PV[,2],
                              SemiPV = ff$fit.BN$ROC$SemiPV[,2]),
              Skewed = data.frame(BN = ff$fit.Skewed$ROC$BN[,2],
                                  PV = ff$fit.Skewed$ROC$PV[,2],
                                  SemiPV = ff$fit.Skewed$ROC$SemiPV[,2]),
              Mixed = data.frame(BN = ff$fit.Mixed$ROC$BN[,2],
                                 PV = ff$fit.Mixed$ROC$PV[,2],
                                 SemiPV = ff$fit.Mixed$ROC$SemiPV[,2])), along = 3)
  ROC1 = abind(list(BN = data.frame(BN = ff$fit.BN$ROC$BN[,3],
                              PV = ff$fit.BN$ROC$PV[,3],
                              SemiPV = ff$fit.BN$ROC$SemiPV[,3]),
              Skewed = data.frame(BN = ff$fit.Skewed$ROC$BN[,3],
                                  PV = ff$fit.Skewed$ROC$PV[,3],
                                  SemiPV = ff$fit.Skewed$ROC$SemiPV[,3]),
              Mixed = data.frame(BN = ff$fit.Mixed$ROC$BN[,3],
                                 PV = ff$fit.Mixed$ROC$PV[,3],
                                 SemiPV = ff$fit.Mixed$ROC$SemiPV[,3])), along = 3)

  J0 = as.matrix(rbind(ff$fit.BN$J.tab$cutoff0[1,],
                       ff$fit.Skewed$J.tab$cutoff0[1,],
                       ff$fit.Mixed$J.tab$cutoff0[1,]))
  J1 = as.matrix(rbind(ff$fit.BN$J.tab$cutoff1[1,],
                       ff$fit.Skewed$J.tab$cutoff1[1,],
                       ff$fit.Mixed$J.tab$cutoff1[1,]))
  rownames(J0) = rownames(J1) = c("BN","Skewed","Mixed")
  
  ER0 = as.matrix(rbind(ff$fit.BN$ER.tab$cutoff0[1,],
                        ff$fit.Skewed$ER.tab$cutoff0[1,],
                        ff$fit.Mixed$ER.tab$cutoff0[1,]))
  ER1 = as.matrix(rbind(ff$fit.BN$ER.tab$cutoff1[1,],
                        ff$fit.Skewed$ER.tab$cutoff1[1,],
                        ff$fit.Mixed$ER.tab$cutoff1[1,]))
  rownames(ER0) = rownames(ER1) = c("BN","Skewed","Mixed")
  
  CZ0 = as.matrix(rbind(ff$fit.BN$CZ.tab$cutoff0[1,],
                        ff$fit.Skewed$CZ.tab$cutoff0[1,],
                        ff$fit.Mixed$CZ.tab$cutoff0[1,]))
  CZ1 = as.matrix(rbind(ff$fit.BN$CZ.tab$cutoff1[1,],
                        ff$fit.Skewed$CZ.tab$cutoff1[1,],
                        ff$fit.Mixed$CZ.tab$cutoff1[1,]))
  rownames(CZ0) = rownames(CZ1) = c("BN","Skewed","Mixed")
  
  IU0 = as.matrix(rbind(ff$fit.BN$IU.tab$cutoff0[1,],
                        ff$fit.Skewed$IU.tab$cutoff0[1,],
                        ff$fit.Mixed$IU.tab$cutoff0[1,]))
  IU1 = as.matrix(rbind(ff$fit.BN$IU.tab$cutoff1[1,],
                        ff$fit.Skewed$IU.tab$cutoff1[1,],
                        ff$fit.Mixed$IU.tab$cutoff1[1,]))
  rownames(IU0) = rownames(IU1) = c("BN","Skewed","Mixed")
  
  AUC0.tab[,,k] = AUC0; AUC1.tab[,,k] = AUC1
  J0.tab[,,k] = J0; J1.tab[,,k] = J1
  ER0.tab[,,k] = ER0; ER1.tab[,,k] = ER1
  CZ0.tab[,,k] = CZ0; CZ1.tab[,,k] = CZ1
  IU0.tab[,,k] = IU0; IU1.tab[,,k] = IU1
  ROC0.tab[,,,k] = ROC0
  ROC1.tab[,,,k] = ROC1

}
Sys.time() - tt

## True values

TrueAUC = rbind(ff$fit.BN$TrueAUC,
                ff$fit.Skewed$TrueAUC,
                ff$fit.Mixed$TrueAUC)
colnames(TrueAUC) = c("x=0", "x=1")
rownames(TrueAUC) = c("BN","Skewed","Mixed")
TrueROC0 = cbind(ff$fit.BN$TrueROC[,1],
                 ff$fit.Skewed$TrueROC[,1],
                 ff$fit.Mixed$TrueROC[,1])
TrueROC1 = cbind(ff$fit.BN$TrueROC[,2],
                 ff$fit.Skewed$TrueROC[,2],
                 ff$fit.Mixed$TrueROC[,2])
colnames(TrueROC0) = colnames(TrueROC1) = c("BN","Skewed","Mixed")
TrueROC = list(`x=0` = TrueROC0, `x=1` = TrueROC1)
TrueCutoff0 = cbind(ff$fit.BN$True.cutoff$cutoff0[,1],
                    ff$fit.Skewed$True.cutoff$cutoff0[,1],
                    ff$fit.Mixed$True.cutoff$cutoff0[,1])
TrueCutoff1 = cbind(ff$fit.BN$True.cutoff$cutoff1[,1],
                    ff$fit.Skewed$True.cutoff$cutoff1[,1],
                    ff$fit.Mixed$True.cutoff$cutoff1[,1])
rownames(TrueCutoff0) = rownames(TrueCutoff1) = c("J","ER","CZ","IU")
colnames(TrueCutoff0) = colnames(TrueCutoff1) = c("BN","Skewed","Mixed")
TrueCutoff = list(`x=0` = TrueCutoff0, `x=1` = TrueCutoff1)

TrueAUC
TrueROC
TrueCutoff

## Estimated values

AUC0.hat = apply(AUC0.tab, 1:2, mean)
ROC0.hat = apply(ROC0.tab, 1:3, mean)
J0.hat = apply(J0.tab, 1:2, mean)
ER0.hat = apply(ER0.tab, 1:2, mean)
CZ0.hat = apply(CZ0.tab, 1:2, mean)
IU0.hat = apply(IU0.tab, 1:2, mean)

AUC1.hat = apply(AUC1.tab, 1:2, mean)
ROC1.hat = apply(ROC1.tab, 1:3, mean)
J1.hat = apply(J1.tab, 1:2, mean)
ER1.hat = apply(ER1.tab, 1:2, mean)
CZ1.hat = apply(CZ1.tab, 1:2, mean)
IU1.hat = apply(IU1.tab, 1:2, mean)

setwd("../")
