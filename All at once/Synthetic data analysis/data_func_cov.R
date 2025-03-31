
# rm(list = ls())
# setwd("~/Library/CloudStorage/OneDrive-UniversityofVirginia/Methodological Research/Optimal Cutoff ROC/Simulation")
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

####################################
###### Function estimating area
###### under the curve numerically
###################################

trap <- function(y, grid = seq(0,1,length=100)){
  # grid <- seq(0,1,length=100)
  val <- trapz(grid,y)
  return(val)
}

####
#### Optimal cut-point function
#### Depends on se and sp
####

J = function(c,x){
  ## negative because maximizing in optim
  -(se(c,x) + sp(c,x) - 1)
}
ER = function(c,x){
  sqrt( (1-se(c))^2 + (1-sp(c))^2 )
}
CZ = function(c,x){
  ## negative because maximizing in optim
  -( se(c,x)*sp(c,x) )
}

IU = function(c,x){
  abs(se(c,x) - mean(AUC)) + abs(sp(c,x) - mean(AUC))
}

####################
####
####  Model fitting
####  function
####
####################


BiNormal_JAGS_fit <- function(y0, y1, X0, X1, burnin = 5000, numb.iter = 5000){
  
  #### JAGS code
  jags_binorm_ROC = "
  model {
  for (i in 1:n0) {
  y0[i] ~ dnorm(b00 + b01*x0[i], tau0)
  }
  for (i in 1:n1) {
  y1[i] ~ dnorm(b10 + b11*x1[i], tau1)
  }
  
  ###################
  ## Priors #########
  ###################
  
  b00 ~ dnorm(0,1/1000)
  b10 ~ dnorm(0,1/1000)
  b01 ~ dnorm(0,1/1000)
  b11 ~ dnorm(0,1/1000)
  tau0 ~ dgamma(0.01,0.01)
  sigma0 <- 1/sqrt(tau0)
  mu1 ~ dnorm(0,1/100)
  tau1 ~ dgamma(0.01,0.01)
  sigma1 <- 1/sqrt(tau1)
  
  }"
  
  N0 <- length(y0)
  N1 <- length(y1)  
  #########################
  #### FIt JAGS model
  #########################
  dat = list(n0=N0, n1=N1, y0=y0, y1=y1, x0 = X0, x1 = X1)
  binorm.model = jags.model(textConnection(jags_binorm_ROC), data = dat, n.chains = 3, quiet = TRUE)
  update(binorm.model, burnin, progress.bar="none")
  
  binorm = jags.samples(binorm.model, variable.names=c("b00","b10","b01","b11","sigma0","sigma1"), 
                        numb.iter,thin=numb.iter/5000, progress.bar="none")

  #########################
  ##### JAGS estimates
  #########################
  b00.hat <- apply(binorm$b00,1:2,mean)
  b10.hat <- apply(binorm$b10,1:2,mean)
  b01.hat <- apply(binorm$b01,1:2,mean)
  b11.hat <- apply(binorm$b11,1:2,mean)
  sigma0.hat <- apply(binorm$sigma0,1:2,mean)
  sigma1.hat <- apply(binorm$sigma1,1:2,mean)
  
  est.tab = list(b00 = b00.hat, b10 = b10.hat,
                       b01 = b01.hat, b11 = b11.hat,
                       sigma0 = sigma0.hat,
                       sigma1 = sigma1.hat)
  return(est.tab)
}

BN_accuracy_cutoff <- function(y0, y1, b00, b01, b10, b11, sigma0, sigma1, x){
  
  J = function(c,x){
    ## negative because maximizing in optim
    -(se(c,x) + sp(c,x) - 1)
  }
  ER = function(c,x){
    sqrt( (1-se(c,x))^2 + (1-sp(c,x))^2 )
  }
  CZ = function(c,x){
    ## negative because maximizing in optim
    -( se(c,x)*sp(c,x) )
  }
  
  IU = function(c,x){
    abs(se(c,x) - mean(AUC)) + abs(sp(c,x) - mean(AUC))
  }
  
  mu0.hat = b00 + b01*x
  mu1.hat = b10 + b11*x
  a.hat = (mu1.hat-mu0.hat)/sigma1
  b.hat = sigma0/sigma1
  grid = seq(0,1,length=100)
  
  ROC <- pnorm(a.hat + b.hat*qnorm(grid))
  AUC <- pnorm(a.hat/sqrt(1+b.hat^2))
  
  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  se = function(c,x){
    1 - pnorm(c, b10+b11*x, sigma1)
  }
  sp = function(c,x){
    pnorm(c, b00+b01*x, sigma0)
  }
  
  J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS", x = x)
  ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS", x = x)
  CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS", x = x)
  IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS", x = x)
  
  J.est = c(J.c$par, -J.c$value, se(J.c$par, x = x), sp(J.c$par, x = x), mean(y1 > J.c$par), mean(y0 < J.c$par))
  ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = x), sp(ER.c$par, x = x), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
  CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = x), sp(CZ.c$par, x = x), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
  IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = x), sp(IU.c$par, x = x), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

BiNormal_reg_fit <- function(y0, y1, X0, X1, x, burnin = 5000, numb.iter = 5000){
  
  fit = BiNormal_JAGS_fit(y0=y0, y1=y1, X0=X0, X1=X1,
                          burnin = 5000, numb.iter = 5000)
  b00 = fit$b00
  b01 = fit$b01
  b10 = fit$b10
  b11 = fit$b11
  sigma0 = fit$sigma0
  sigma1 = fit$sigma1
  
  grid = seq(0, 1, length = 100)
  G = length(b00)
  # x = c(60,80)
  
  AUC.x0 = AUC.x1 = rep(NA, G)
  ROC.x0 = ROC.x1 = array(NA, dim = c(length(grid), G))
  cutoff.x0 = cutoff.x1 = array(NA, dim = c(4, G))
  
  for(g in 1:G){
    bn0 = BN_accuracy_cutoff(y0, y1, b00=b00[g], b01=b01[g], b10=b10[g], b11=b11[g], 
                             sigma0=sigma0[g], sigma1=sigma1[g], x=x[1])
    bn1 = BN_accuracy_cutoff(y0, y1, b00=b00[g], b01=b01[g], b10=b10[g], b11=b11[g], 
                             sigma0=sigma0[g], sigma1=sigma1[g], x=x[2])
    AUC.x0[g] = bn0$AUC
    AUC.x1[g] = bn1$AUC
    ROC.x0[,g] = bn0$ROC
    ROC.x1[,g] = bn1$ROC
    cutoff.x0[,g] = bn0$cutoff.tab[,1]
    cutoff.x1[,g] = bn1$cutoff.tab[,1]
  }
  
  ROC0.hat = apply(ROC.x0, 1, mean)
  ROC1.hat = apply(ROC.x1, 1, mean)
  AUC0.hat = c(mean(AUC.x0), sd(AUC.x0), quantile(AUC.x0, probs = c(0.025, 0.975)))
  AUC0.hat = matrix(AUC0.hat, nrow = 1, byrow = TRUE)
  AUC1.hat = c(mean(AUC.x1), sd(AUC.x1), quantile(AUC.x1, probs = c(0.025, 0.975)))
  AUC1.hat = matrix(AUC1.hat, nrow = 1, byrow = TRUE)
  colnames(AUC0.hat) = colnames(AUC1.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  cutoff0.hat = cbind(apply(cutoff.x0, 1, mean),
                     apply(cutoff.x0, 1, sd),
                     t(apply(cutoff.x0, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  cutoff1.hat = cbind(apply(cutoff.x1, 1, mean),
                      apply(cutoff.x1, 1, sd),
                      t(apply(cutoff.x1, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  rownames(cutoff0.hat) = rownames(cutoff1.hat) = c("J","ER","CZ","IU")
  colnames(cutoff0.hat) = colnames(cutoff1.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  AUC = list(AUC0 = AUC0.hat, AUC1 = AUC1.hat)
  ROC = data.frame(ROC0 = ROC0.hat, ROC1 = ROC1.hat)
  cutoff = list(cutoff0 = cutoff0.hat,
                cutoff1 = cutoff1.hat)
  out = list(AUC = AUC, ROC = ROC, cutoff = cutoff)
  return(out)
}


BayesPV_reg_fit <- function(y0, y1, X0, X1, x, burnin = 5000, numb.iter = 5000){
  
  J = function(c,x){
    ## negative because maximizing in optim
    -(se(c,x) + sp(c,x) - 1)
  }
  ER = function(c,x){
    sqrt( (1-se(c,x))^2 + (1-sp(c,x))^2 )
  }
  CZ = function(c,x){
    ## negative because maximizing in optim
    -( se(c,x)*sp(c,x) )
  }
  
  IU = function(c,x){
    abs(se(c,x) - mean(AUC)) + abs(sp(c,x) - mean(AUC))
  }
  
  healthy_fit_mod = "
  model {
  #model
  for (i in 1:length(y0)){
  y0[i] ~ dnorm(b00 + b01*x0[i],tau0)
  }
  #priors
  b00  ~ dnorm(0, 0.0001)
  b01  ~ dnorm(0, 0.0001)
  tau0 ~ dgamma(0.001,0.001)
  sigma0 <- 1/sqrt(tau0)
  }"
  
  PV_fit_mod = "
  model {
  #model
  for (i in 1:length(z.star)){
  z.star[i] ~ dnorm(b0 + b1*x1[i],tau)
  }
  #priors
  b0  ~ dnorm(0, 0.0001)
  b1  ~ dnorm(0, 0.0001)
  tau ~ dgamma(0.001,0.001)
  sigma <- 1/sqrt(tau)
  }"
  
  grid <- seq(0,1,length=100)
  
  ##stage 1
  
  tm <- Sys.time()
  dat_stage1 = list(y0 = y0, x0 = X0)
  stage1.model = jags.model(textConnection(healthy_fit_mod), data = dat_stage1, 
                            n.chains = 3, quiet = TRUE)
  update(stage1.model, burnin, progress.bar = "none")
  stage1_jags = jags.samples(stage1.model, variable.names=c("b00", "b01", "sigma0"), 
                             numb.iter,thin=numb.iter/5000, progress.bar = "none")
  Sys.time()-tm
  
  b00 <- as.vector(apply(stage1_jags$b00,1:2,mean))
  b01 <- as.vector(apply(stage1_jags$b01,1:2,mean))
  sigma0 <- as.vector(apply(stage1_jags$sigma0,1:2,mean))
  b00.hat = mean(b00)
  b01.hat = mean(b01)
  sigma0.hat = mean(sigma0)
  
  PV <- 1-pnorm(y1, b00.hat + b01.hat*X1, sigma0.hat)
  PV <- ifelse(PV<(1e-15),(1e-15),ifelse(PV>(1-(1e-15)),(1-(1e-15)),PV))
  PV_one = PV
  
  #####################################
  #####################################
  #######                       #######
  #######    PV Model           #######
  #######                       #######
  #####################################
  #####################################
  
  tm <- Sys.time()
  dat_normal_reg = list(z.star=qnorm(PV_one), x1 = X1)
  normal.model = jags.model(textConnection(PV_fit_mod), 
                            data = dat_normal_reg, n.chains = 3, quiet = TRUE)
  update(normal.model, burnin, progress.bar = "none")
  normal_reg_jags = jags.samples(normal.model, variable.names=c("b0","b1","sigma"), 
                                 numb.iter,thin=numb.iter/5000, progress.bar = "none")
  Sys.time()-tm
  
  b0 <- as.vector(apply(normal_reg_jags$b0,1:2,mean)) 
  b1 <- as.vector(apply(normal_reg_jags$b1,1:2,mean)) 
  sigma <- as.vector(apply(normal_reg_jags$sigma,1:2,mean))  

  G = length(b0)
  # x = c(60,80)
  
  ROC.x0 = sapply(1:G, function(g){pnorm(qnorm(grid), b0[g] + b1[g]*x[1], sigma[g])})
  ROC.x1 = sapply(1:G, function(g){pnorm(qnorm(grid), b0[g] + b1[g]*x[2], sigma[g])})
  AUC.x0 = apply(ROC.x0,2,trap)
  AUC.x1 = apply(ROC.x1,2,trap)

  ROC0.hat = apply(ROC.x0, 1, mean)
  ROC1.hat = apply(ROC.x1, 1, mean)
  AUC0.hat = c(mean(AUC.x0), sd(AUC.x0), quantile(AUC.x0, probs = c(0.025, 0.975)))
  AUC0.hat = matrix(AUC0.hat, nrow = 1, byrow = TRUE)
  AUC1.hat = c(mean(AUC.x1), sd(AUC.x1), quantile(AUC.x1, probs = c(0.025, 0.975)))
  AUC1.hat = matrix(AUC1.hat, nrow = 1, byrow = TRUE)
  colnames(AUC0.hat) = colnames(AUC1.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  
  PV0.est.func = function(g){
    
    ROC = pnorm(qnorm(grid), b0[g] + b1[g]*x[1], sigma[g])
    AUC = trap(ROC)

    ##########################
    ###########
    ##### Optimal functions
    ###########
    ##########################
    
    J = function(c){
      ## negative because maximizing in optim
      -(se(c) + sp(c) - 1)
    }
    ER = function(c){
      sqrt( (1-se(c))^2 + (1-sp(c))^2 )
    }
    CZ = function(c){
      ## negative because maximizing in optim
      -( se(c)*sp(c) )
    }
    
    IU = function(c){
      abs(se(c) - mean(AUC)) + abs(sp(c) - mean(AUC))
    }
    
    se = function(c){
      pnorm(qnorm(1 - pnorm(c, b00[g] + b01[g]*x[1], sigma0[g])), b0[g] + b1[g]*x[1], sigma[g])
    }
    
    sp = function(c){
      pnorm(c, b00[g] + b01[g]*x[1], sigma0[g])
    }
    
    J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS")
    ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS")
    CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS")
    IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS")
    
    J.est = c(J.c$par)
    ER.est = c(ER.c$par)
    CZ.est = c(CZ.c$par)
    IU.est = c(IU.c$par)
    cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
    return(cutoff = cutoff.tab)
  }
  PV1.est.func = function(g){
    
    ROC = pnorm(qnorm(grid), b0[g] + b1[g]*x[2], sigma[g])
    AUC = trap(ROC)
    
    ##########################
    ###########
    ##### Optimal functions
    ###########
    ##########################
    
    J = function(c){
      ## negative because maximizing in optim
      -(se(c) + sp(c) - 1)
    }
    ER = function(c){
      sqrt( (1-se(c))^2 + (1-sp(c))^2 )
    }
    CZ = function(c){
      ## negative because maximizing in optim
      -( se(c)*sp(c) )
    }
    
    IU = function(c){
      abs(se(c) - mean(AUC)) + abs(sp(c) - mean(AUC))
    }
    
    se = function(c){
      pnorm(qnorm(1 - pnorm(c, b00[g] + b01[g]*x[2], sigma0[g])), b0[g] + b1[g]*x[2], sigma[g])
    }
    
    sp = function(c){
      pnorm(c, b00[g] + b01[g]*x[2], sigma0[g])
    }
    
    J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS")
    ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS")
    CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS")
    IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS")
    
    J.est = c(J.c$par)
    ER.est = c(ER.c$par)
    CZ.est = c(CZ.c$par)
    IU.est = c(IU.c$par)
    cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
    return(cutoff = cutoff.tab)
  }

  cutoff0.all = sapply(1:G, PV0.est.func)
  cutoff1.all = sapply(1:G, PV1.est.func)
  cutoff0.hat = cbind(apply(cutoff0.all, 1, mean),
                     apply(cutoff0.all, 1, sd),
                     t(apply(cutoff0.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  cutoff1.hat = cbind(apply(cutoff1.all, 1, mean),
                      apply(cutoff1.all, 1, sd),
                      t(apply(cutoff1.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  
  rownames(cutoff0.hat) = rownames(cutoff1.hat) = c("J","ER","CZ","IU")
  colnames(cutoff0.hat) = colnames(cutoff1.hat) = c("Mean","SD","Low 95%", "Upper 95%")

  AUC = list(AUC0 = AUC0.hat, AUC1 = AUC1.hat)
  ROC = data.frame(ROC0 = ROC0.hat, ROC1 = ROC1.hat)
  cutoff = list(cutoff0 = cutoff0.hat,
                cutoff1 = cutoff1.hat)
  out = list(AUC = AUC, ROC = ROC, cutoff = cutoff)
  return(out)
}

SemiPar_reg_fit <- function(y0, y1, X0, X1, x, burnin = 5000, numb.iter = 5000){
  
  J = function(c,x){
    ## negative because maximizing in optim
    -(se(c,x) + sp(c,x) - 1)
  }
  ER = function(c,x){
    sqrt( (1-se(c,x))^2 + (1-sp(c,x))^2 )
  }
  CZ = function(c,x){
    ## negative because maximizing in optim
    -( se(c,x)*sp(c,x) )
  }
  
  IU = function(c,x){
    abs(se(c,x) - mean(AUC)) + abs(sp(c,x) - mean(AUC))
  }
  
  dp_stage1 = "
    model {
    for (i in 1:n0) {
    y0[i] ~ dnorm(b00[zeta0[i]] + b01[zeta0[i]]*x0[i], tau0)
    zeta0[i] ~ dcat(pi0[])
    }
    
    ###################
    ## Priors #########
    ###################
    for (h in 1:H) {
    b00[h] ~ dnorm(0, 0.0001)
    b01[h] ~ dnorm(0, 0.0001)
    }
    tau0 ~ dgamma(0.1,0.1)
    sigma0 <- 1/sqrt(tau0)
    
    ########################################
    ## Concentration parameter #############
    ########################################
    a0 ~ dgamma(1,1)

    ##################################
    ## Stick breaking ################
    ##################################
    for (h in 1:(H-1)) { 
    V0[h] ~ dbeta(1,a0)T(1.0E-7,1-1.0E-7) 
    }
    V0[H] <- 1
    pi0[1] <- V0[1]

    for (h in 2:H) {
    pi0[h] <- V0[h] * (1-V0[h-1]) * pi0[h-1] / V0[h-1]
    }
    }"
  
  dp_PV = "
    model {
    for (i in 1:n1) {
    z[i] ~ dnorm(b0[zeta[i]] + b1[zeta[i]]*x1[i], tau)
    #logit(eta[i]) <- z[i]
    }
    for (i in 1:n1) {
    zeta[i] ~ dcat(pi[])
    }
    for (h in 1:H) {
    b0[h] ~ dnorm(0, 0.0001)
    b1[h] ~ dnorm(0, 0.0001)
    }
    tau ~ dgamma(0.1,0.1)
    sigma <- 1/sqrt(tau)

    a ~ dgamma(1,1)
    
    # Stick breaking
    for (h in 1:(H-1)) { 
    V[h] ~ dbeta(1,a)T(1.0E-7,1-1.0E-7) 
    }
    V[H] <- 1
    pi[1] <- V[1]
    for (h in 2:H) {
    pi[h] <- V[h] * (1-V[h-1]) * pi[h-1] / V[h-1]
    }
    }"
  
  #########################
  #### Stage 1
  #########################
  
  #########################
  #### FIt JAGS model
  #########################
  dat = list(n0 = length(y0), H=30, y0 = y0, x0 = X0)
  CDF.model = jags.model(textConnection(dp_stage1), data = dat, n.chains = 3, quiet=TRUE)
  update(CDF.model, burnin, progress.bar="none")
  CDF = jags.samples(CDF.model, variable.names=c("b00","b01","sigma0","pi0"), 
                     numb.iter,thin=numb.iter/5000, progress.bar="none")

  #########################
  ##### JAGS estimates
  #########################
  b00.hat <- apply(CDF$b00,1:2,mean)
  b01.hat <- apply(CDF$b01,1:2,mean)
  sigma0.hat <- apply(CDF$sigma0,1:2,mean)
  pi0.hat <- apply(CDF$pi0,1:2,mean)
  G <- dim(b00.hat)[2]
  
  ########################
  ###### Estimate PV
  ########################
  
  n1 = length(y1)
  PV.tab <-  rep(NA, n1)
  for(j in 1:n1)
  {
    # PV.tab[j] <- 1-sum(pi0.hat*pnorm(y1[j], b00.hat + b01.hat*X1[j], sigma0.hat))
    PV.tab[j] <- 1-sum(apply(pi0.hat,1,mean)*pnorm(y1[j], apply(b00.hat,1,mean) + apply(b01.hat,1,mean)*X1[j], apply(sigma0.hat,1,mean)))
  }
  eta = PV.tab
  
  eta = ifelse(eta < 0.000000001, 0.000000001, 
               ifelse(eta > (1-0.000000001), (1-0.000000001), eta))
  
  #########################
  #### Stage 2
  #########################
  
  #########################
  #### FIt JAGS model
  #########################
  
  dat = list(n1=length(y1), H=30, z=qnorm(eta), x1 = X1)
  PV.model = jags.model(textConnection(dp_PV), data = dat, n.chains = 3, quiet=TRUE)
  update(PV.model, burnin,progress.bar = "none")
  PVjags = jags.samples(PV.model, variable.names=c("b0","b1","sigma","pi","a"), 
                        numb.iter,thin=numb.iter/5000,progress.bar = "none")
  
  #########################
  ##### JAGS estimates
  #########################

  b0.hat <- apply(PVjags$b0,1:2,mean)
  b1.hat <- apply(PVjags$b1,1:2,mean)
  sigma.hat <- apply(PVjags$sigma,1:2,mean)
  pi.hat <- apply(PVjags$pi,1:2,mean)
  G <- dim(b0.hat)[2]
  
  ##################################
  ##### Estimate CDF (ROC) at 
  ##### different iteration level
  ##################################
  grid <- seq(0,1,length=100)

  ROC.x0 = sapply(1:G, function(g){
    
    ROC.tab <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC.tab[j] <- sum(pi.hat[,g]*pnorm(qnorm(grid[j]), b0.hat[,g] + b1.hat[,g]*x[1], sigma.hat[1,g]))
    }
    return(ROC.tab)
    })
  ROC.x1 = sapply(1:G, function(g){
    
    ROC.tab <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC.tab[j] <- sum(pi.hat[,g]*pnorm(qnorm(grid[j]), b0.hat[,g] + b1.hat[,g]*x[2], sigma.hat[1,g]))
    }
    return(ROC.tab)
  })
  AUC.x0 = apply(ROC.x0,2,trap)
  AUC.x1 = apply(ROC.x1,2,trap)
  
  ROC0.hat = apply(ROC.x0, 1, mean)
  ROC1.hat = apply(ROC.x1, 1, mean)
  AUC0.hat = c(mean(AUC.x0), sd(AUC.x0), quantile(AUC.x0, probs = c(0.025, 0.975)))
  AUC0.hat = matrix(AUC0.hat, nrow = 1, byrow = TRUE)
  AUC1.hat = c(mean(AUC.x1), sd(AUC.x1), quantile(AUC.x1, probs = c(0.025, 0.975)))
  AUC1.hat = matrix(AUC1.hat, nrow = 1, byrow = TRUE)
  colnames(AUC0.hat) = colnames(AUC1.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  SemiPV0.est.func = function(g){
    
      ROC <-  rep(NA, length(grid))
      for(j in 1:length(grid))
      {
        ROC[j] <- sum(pi.hat[,g]*pnorm(qnorm(grid[j]), b0.hat[,g] + b1.hat[,g]*x[1], sigma.hat[1,g]))
      }
      AUC = trap(ROC)
      
      ##########################
      ###########
      ##### Optimal functions
      ###########
      ##########################
      
      J = function(c){
        ## negative because maximizing in optim
        -(se(c) + sp(c) - 1)
      }
      ER = function(c){
        sqrt( (1-se(c))^2 + (1-sp(c))^2 )
      }
      CZ = function(c){
        ## negative because maximizing in optim
        -( se(c)*sp(c) )
      }
      
      IU = function(c){
        abs(se(c) - mean(AUC)) + abs(sp(c) - mean(AUC))
      }
      
      sp = function(c){
        sum(pi0.hat[,g]*pnorm(c, b00.hat[,g] + b01.hat[,g]*x[1], sigma0.hat[1,g]))
      }
      
      se = function(c){
        sum(pi.hat[,g]*pnorm(qnorm(1-sum(pi0.hat[,g]*pnorm(c, b00.hat[,g] + b01.hat[,g]*x[1], sigma0.hat[1,g]))), 
                             b0.hat[,g] + b1.hat[,g]*x[1], sigma.hat[1,g]))
      }
      
      J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS")
      ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS")
      CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS")
      IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS")
      
      J.est = c(J.c$par)
      ER.est = c(ER.c$par)
      CZ.est = c(CZ.c$par)
      IU.est = c(IU.c$par)
      cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
      return(cutoff = cutoff.tab)
  }
  SemiPV1.est.func = function(g){
    
    ROC <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC[j] <- sum(pi.hat[,g]*pnorm(qnorm(grid[j]), b0.hat[,g] + b1.hat[,g]*x[2], sigma.hat[1,g]))
    }
    AUC = trap(ROC)
    
    ##########################
    ###########
    ##### Optimal functions
    ###########
    ##########################
    
    J = function(c){
      ## negative because maximizing in optim
      -(se(c) + sp(c) - 1)
    }
    ER = function(c){
      sqrt( (1-se(c))^2 + (1-sp(c))^2 )
    }
    CZ = function(c){
      ## negative because maximizing in optim
      -( se(c)*sp(c) )
    }
    
    IU = function(c){
      abs(se(c) - mean(AUC)) + abs(sp(c) - mean(AUC))
    }
    
    sp = function(c){
      sum(pi0.hat[,g]*pnorm(c, b00.hat[,g] + b01.hat[,g]*x[2], sigma0.hat[1,g]))
    }
    
    se = function(c){
      sum(pi.hat[,g]*pnorm(qnorm(1-sum(pi0.hat[,g]*pnorm(c, b00.hat[,g] + b01.hat[,g]*x[2], sigma0.hat[1,g]))), 
                           b0.hat[,g] + b1.hat[,g]*x[2], sigma.hat[1,g]))
    }
    
    J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS")
    ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS")
    CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS")
    IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS")
    
    J.est = c(J.c$par)
    ER.est = c(ER.c$par)
    CZ.est = c(CZ.c$par)
    IU.est = c(IU.c$par)
    cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
    return(cutoff = cutoff.tab)
  }
  
  cutoff0.all = sapply(1:G, SemiPV0.est.func)
  cutoff1.all = sapply(1:G, SemiPV1.est.func)
  cutoff0.hat = cbind(apply(cutoff0.all, 1, mean),
                      apply(cutoff0.all, 1, sd),
                      t(apply(cutoff0.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  cutoff1.hat = cbind(apply(cutoff1.all, 1, mean),
                      apply(cutoff1.all, 1, sd),
                      t(apply(cutoff1.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  
  rownames(cutoff0.hat) = rownames(cutoff1.hat) = c("J","ER","CZ","IU")
  colnames(cutoff0.hat) = colnames(cutoff1.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  AUC = list(AUC0 = AUC0.hat, AUC1 = AUC1.hat)
  ROC = data.frame(ROC0 = ROC0.hat, ROC1 = ROC1.hat)
  cutoff = list(cutoff0 = cutoff0.hat,
                cutoff1 = cutoff1.hat)
  out = list(AUC = AUC, ROC = ROC, cutoff = cutoff)
  return(out)
}

