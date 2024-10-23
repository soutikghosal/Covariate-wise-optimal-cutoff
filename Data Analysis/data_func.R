

# rm(list=ls())
# setwd("~/Library/CloudStorage/OneDrive-UniversityofVirginia/Methodological Research/Optimal Cutoff ROC/Simulation")
library("rjags")
library(edgeR)
library(dplyr)
library(plyr)
library(ggplot2)
library(coda)
require(devtools)
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

####################
####
####  Model fitting
####  function
####
####################

Empirical_func <- function(y0, y1){
  
  n0 = length(y0)
  n1 = length(y1)
  
  F0 = function(c){
    # (1/n0)*sum((y0 <= c))
    ecdf(y0)(c)
  }
  F0.inv = function(c){
    quantile(y0, c)
  }
  F1 = function(c){
    # (1/n1)*sum((y1 <= c))
    ecdf(y1)(c)
  }
  
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
  
  grid = seq(0,1,length=100)
  ROC = 1-F1(F0.inv(1-grid))
  AUC = trap(ROC)
  
  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  se = function(c){
    1 - F1(c)
  }
  sp = function(c){
    F0(c)
  }
  
  J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
  ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
  CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
  IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
  
  J.est = c(J.c$par)
  ER.est = c(ER.c$par)
  CZ.est = c(CZ.c$par)
  IU.est = c(IU.c$par)
  cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
  
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

## Bootstrap version of empirical to calculate SE
Empirical_Boot <- function(y0, y1, B = 1000){
  
  # df = data.frame(y = c(y0,y1), d = c(rep(0, length(y0)),rep(1, length(y1))))
  df0 = data.frame(y0)
  df1 = data.frame(y1)
  
  grid = seq(0,1,length=100)
  ROC.tab = array(NA, dim = c(length(grid), B))
  AUC.tab = rep(NA, B)
  cutoff.all = array(NA, dim = c(4, B))
  
  for(b in 1:B){
    id.y0 = sample.int(length(y0), size = length(y0), replace = TRUE)
    id.y1 = sample.int(length(y1), size = length(y1), replace = TRUE)
    
    y0.new = df0[id.y0,]
    y1.new = df1[id.y1,]
    
    emp.fit = Empirical_func(y0.new, y1.new)
    AUC.tab[b] = emp.fit$AUC
    ROC.tab[,b] = emp.fit$ROC
    cutoff.all[,b] = emp.fit$cutoff.tab
  }
  
  ROC.hat = apply(ROC.tab, 1, mean)
  AUC.hat = c(mean(AUC.tab), sd(AUC.tab), quantile(AUC.tab, probs = c(0.025, 0.975)))
  AUC.hat = matrix(AUC.hat, nrow = 1, byrow = TRUE)
  colnames(AUC.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  # cutoff.all = -cutoff.all
  cutoff.hat = cbind(apply(cutoff.all, 1, mean),
                     apply(cutoff.all, 1, sd),
                     t(apply(cutoff.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  rownames(cutoff.hat) = c("J","ER","CZ","IU")
  colnames(cutoff.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  out=list(ROC=ROC.hat, AUC=AUC.hat, cutoff=cutoff.hat)
  return(out)
}

NonPar_func <- function(y0, y1){
  
  n0 = length(y0)
  n1 = length(y1)
  h0 = 0.9*min(sd(y0), abs(diff(as.numeric(quantile(y0, probs = c(0.75, 0.25)))))/1.34 )*(n0^(-1.5))
  h1 = 0.9*min(sd(y1), abs(diff(as.numeric(quantile(y1, probs = c(0.75, 0.25)))))/1.34 )*(n1^(-1.5))
  
  F0 = function(c){
    CDF = (1/n0)*sum(pnorm((c-y0)/h0))
    return(CDF)
  }
  
  inverse = function(fn, interval = NULL, lower = min(interval), upper = max(interval), ...){
    Vectorize(function(y){
      uniroot(f=function(x){fn(x)-y}, lower=lower, upper=upper, ...)$root
    })
  }
  
  F0.inv = inverse(F0, lower=-100, upper=100)
  F1 = function(c){
    CDF = (1/n1)*sum(pnorm((c-y1)/h1))
    return(CDF)
  }
  
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
  
  grid = seq(0,1,length=100)
  # ROC = 1-F1(F0.inv(1-grid))
  ROC = 1- sapply(F0.inv(1-grid),F1)
  # AUC = trap(ROC)
  AUC = (sum(pnorm(c(outer(y1, y0, "-"))/sqrt(h0^2 + h1^2))))/(n0*n1)
  
  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  se = function(c){
    1 - F1(c)
  }
  sp = function(c){
    F0(c)
  }
  
  J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
  ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
  CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
  IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
  
  J.est = c(J.c$par)
  ER.est = c(ER.c$par)
  CZ.est = c(CZ.c$par)
  IU.est = c(IU.c$par)
  cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
  # out = list(AUC = AUC, ROC = ROC, cutoff = cutoff.tab)
  # return(out)
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

## Bootstrap version of NonPar to calculate SE
NonPar_Boot <- function(y0, y1, B = 1000){
  
  # df = data.frame(y = c(y0,y1), d = c(rep(0, length(y0)),rep(1, length(y1))))
  df0 = data.frame(y0)
  df1 = data.frame(y1)
  
  grid = seq(0,1,length=100)
  ROC.tab = array(NA, dim = c(length(grid), B))
  AUC.tab = rep(NA, B)
  cutoff.all = array(NA, dim = c(4, B))
  
  for(b in 1:B){
    id.y0 = sample.int(length(y0), size = length(y0), replace = TRUE)
    id.y1 = sample.int(length(y1), size = length(y1), replace = TRUE)
    
    y0.new = df0[id.y0,]
    y1.new = df1[id.y1,]
    
    nonpar.fit = NonPar_func(y0.new, y1.new)
    AUC.tab[b] = nonpar.fit$AUC
    ROC.tab[,b] = nonpar.fit$ROC
    cutoff.all[,b] = nonpar.fit$cutoff.tab
  }
  
  ROC.hat = apply(ROC.tab, 1, mean)
  AUC.hat = c(mean(AUC.tab), sd(AUC.tab), quantile(AUC.tab, probs = c(0.025, 0.975)))
  AUC.hat = matrix(AUC.hat, nrow = 1, byrow = TRUE)
  colnames(AUC.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  # cutoff.all = -cutoff.all
  cutoff.hat = cbind(apply(cutoff.all, 1, mean),
                     apply(cutoff.all, 1, sd),
                     t(apply(cutoff.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  rownames(cutoff.hat) = c("J","ER","CZ","IU")
  colnames(cutoff.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  out=list(ROC=ROC.hat, AUC=AUC.hat, cutoff=cutoff.hat)
  return(out)
}

BiNormal_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  #### JAGS code
  jags_binorm_ROC = "
  model {
  for (i in 1:n0) {
  y0[i] ~ dnorm(mu0, tau0)
  }
  for (i in 1:n1) {
  y1[i] ~ dnorm(mu1, tau1)
  }
  
  ###################
  ## Priors #########
  ###################
  
  mu0 ~ dnorm(0,1/100)
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
  dat = list(n0=N0, n1=N1, y0=y0, y1=y1)
  binorm.model = jags.model(textConnection(jags_binorm_ROC), data = dat, n.chains = 3, quiet = TRUE)
  update(binorm.model, burnin, progress.bar="none")
  
  binorm = jags.samples(binorm.model, variable.names=c("mu0","mu1","sigma0","sigma1"), 
                        numb.iter,thin=numb.iter/5000, progress.bar="none")

  #########################
  ##### JAGS estimates
  #########################
  mu0.hat <- apply(binorm$mu0,1:2,mean)
  mu1.hat <- apply(binorm$mu1,1:2,mean)
  sigma0.hat <- apply(binorm$sigma0,1:2,mean)
  sigma1.hat <- apply(binorm$sigma1,1:2,mean)
  
  G = length(mu0.hat)

  ab.est = function(x){
    a.hat = (mu1.hat[x]-mu0.hat[x])/sigma1.hat[x]
    b.hat = sigma0.hat[x]/sigma1.hat[x]
    return(c(a.hat, b.hat))
  }
  BN.ROC.est = function(x){
    grid <- seq(0,1,length=100)
    a.hat = ab.est(x)[1]
    b.hat = ab.est(x)[2]
    ROC <- pnorm(a.hat + b.hat*qnorm(grid))
    return(ROC)
  }
  BN.AUC.est = function(x){
    a.hat = ab.est(x)[1]
    b.hat = ab.est(x)[2]
    AUC <- pnorm(a.hat/sqrt(1+b.hat^2))
    return(AUC)
  }
  ROC.tab = sapply(1:G, BN.ROC.est)
  ROC.hat = apply(ROC.tab, 1, mean)
  AUC.tab = sapply(1:G, BN.AUC.est)
  AUC.hat = c(mean(AUC.tab), sd(AUC.tab), quantile(AUC.tab, probs = c(0.025, 0.975)))
  AUC.hat = matrix(AUC.hat, nrow = 1, byrow = TRUE)
  colnames(AUC.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  BN.est.func = function(x){
    
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
    
    a.hat = (mu1.hat[x]-mu0.hat[x])/sigma1.hat[x]
    b.hat = sigma0.hat[x]/sigma1.hat[x]
    grid = seq(0,1,length=100)
    ROC <- pnorm(a.hat + b.hat*qnorm(grid))
    AUC <- pnorm(a.hat/sqrt(1+b.hat^2))
    
    ##########################
    ###########
    ##### Optimal functions
    ###########
    ##########################
    
    se = function(c){
      1 - pnorm(c, mu1.hat[x], sigma1.hat[x])
    }
    sp = function(c){
      pnorm(c, mu0.hat[x], sigma0.hat[x])
    }
    
    J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
    ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
    CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
    IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
    
    J.est = c(J.c$par)
    ER.est = c(ER.c$par)
    CZ.est = c(CZ.c$par)
    IU.est = c(IU.c$par)
    cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
    # out = list(AUC = AUC, ROC = ROC, cutoff = cutoff.tab)
    # return(out)
    return(cutoff = cutoff.tab)
  }
  cutoff.all = sapply(1:G, BN.est.func)
  # cutoff.all = -cutoff.all
  cutoff.hat = cbind(apply(cutoff.all, 1, mean),
                     apply(cutoff.all, 1, sd),
                     t(apply(cutoff.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  rownames(cutoff.hat) = c("J","ER","CZ","IU")
  colnames(cutoff.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  out=list(ROC=ROC.hat, AUC=AUC.hat, cutoff=cutoff.hat)
  return(out)
}

PV_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  healthy_fit_mod = "
  model {
  #model
  for (i in 1:length(y0)){
  y0[i] ~ dnorm(mu0,tau0)
  }
  #priors
  mu0  ~ dnorm(0, 0.0001)
  tau0 ~ dgamma(0.001,0.001)
  sigma0 <- 1/sqrt(tau0)
  }"
  
  PV_fit_mod = "
  model {
  #model
  for (i in 1:length(z.star)){
  z.star[i] ~ dnorm(mu,tau)
  }
  #priors
  mu  ~ dnorm(0, 0.0001)
  tau ~ dgamma(0.001,0.001)
  sigma <- 1/sqrt(tau)
  }"
  
  grid <- seq(0,1,length=100)
  
  ##stage 1
  
  tm <- Sys.time()
  dat_stage1 = list(y0 = y0)
  stage1.model = jags.model(textConnection(healthy_fit_mod), data = dat_stage1, 
                            n.chains = 3, quiet = TRUE)
  update(stage1.model, burnin, progress.bar = "none")
  stage1_jags = jags.samples(stage1.model, variable.names=c("mu0", "sigma0"), 
                             numb.iter,thin=numb.iter/5000, progress.bar = "none")
  Sys.time()-tm
  
  mu0 <- as.vector(apply(stage1_jags$mu0,1:2,mean))
  sigma0 <- as.vector(apply(stage1_jags$sigma0,1:2,mean))
  # mu0.hat = mean(mu0)
  # sigma0.hat = mean(sigma0)
  
  PV <- 1-pnorm(y1, mean(mu0), mean(sigma0))
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
  dat_normal_reg = list(z.star=qnorm(PV_one))
  normal.model = jags.model(textConnection(PV_fit_mod), 
                                data = dat_normal_reg, n.chains = 3, quiet = TRUE)
  update(normal.model, burnin, progress.bar = "none")
  normal_reg_jags = jags.samples(normal.model, variable.names=c("mu","sigma"), 
                                 numb.iter,thin=numb.iter/5000, progress.bar = "none")
  Sys.time()-tm
  
  mu <- as.vector(apply(normal_reg_jags$mu,1:2,mean)) 
  sigma <- as.vector(apply(normal_reg_jags$sigma,1:2,mean))  
  # mu = mean(mu)
  # sigma = mean(sigma)
  
  G = length(mu0)
  ROC.tab = sapply(1:G, function(x){pnorm(qnorm(grid),mu[x],sigma[x])})
  ROC.hat = apply(ROC.tab, 1, mean)
  AUC.tab = apply(ROC.tab, 2, trap)
  AUC.hat = c(mean(AUC.tab), sd(AUC.tab), quantile(AUC.tab, probs = c(0.025, 0.975)))
  AUC.hat = matrix(AUC.hat, nrow = 1, byrow = TRUE)
  colnames(AUC.hat) = c("Mean","SD","Low 95%", "Upper 95%")

  # ROC <- pnorm(qnorm(grid),mu,sigma)
  # AUC = trap(ROC)

  PV.est.func = function(x){
    
    ROC = pnorm(qnorm(grid),mu[x],sigma[x])
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
      pnorm(c, mu0[x], sigma0[x])
    }
    se = function(c){
      pnorm(qnorm(1-pnorm(c, mu0[x], sigma0[x])), mu[x], sigma[x])
    }
    
    J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
    ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
    CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
    IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
    
    J.est = c(J.c$par)
    ER.est = c(ER.c$par)
    CZ.est = c(CZ.c$par)
    IU.est = c(IU.c$par)
    cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
    # out = list(AUC = AUC, ROC = ROC, cutoff = cutoff.tab)
    # return(out)
    return(cutoff = cutoff.tab)
  }
  
  cutoff.all = sapply(1:G, PV.est.func)
  # cutoff.all = -cutoff.all
  cutoff.hat = cbind(apply(cutoff.all, 1, mean),
                     apply(cutoff.all, 1, sd),
                     t(apply(cutoff.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  rownames(cutoff.hat) = c("J","ER","CZ","IU")
  colnames(cutoff.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  out=list(ROC=ROC.hat, AUC=AUC.hat, cutoff=cutoff.hat)
  return(out)
}

SemiPar_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  dp_stage1_old = "
    model {
    for (i in 1:n0) {
    y0[i] ~ dnorm(mu0[zeta0[i]], tau0[zeta0[i]])
    zeta0[i] ~ dcat(pi0[])
    }
    
    ###################
    ## Priors #########
    ###################
    for (h in 1:H) {
    mu0[h] ~ dnorm(0, 0.001)
    tau0[h] ~ dgamma(0.1,0.1)
    sigma0[h] <- 1/sqrt(tau0[h])
    }
    
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
  
  dp_stage1 = "
    model {
    for (i in 1:n0) {
    y0[i] ~ dnorm(mu0[zeta0[i]], tau0)
    zeta0[i] ~ dcat(pi0[])
    }
    
    ###################
    ## Priors #########
    ###################
    for (h in 1:H) {
    mu0[h] ~ dnorm(0, 0.001)
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
    z[i] ~ dnorm(mu[zeta[i]], tau)
    #logit(eta[i]) <- z[i]
    }
    for (i in 1:n1) {
    zeta[i] ~ dcat(pi[])
    }
    for (h in 1:H) {
    mu[h] ~ dnorm(0, 0.001)
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
  dat = list(n0 = length(y0), H=30, y0 = y0)
  CDF.model = jags.model(textConnection(dp_stage1), data = dat, n.chains = 3, quiet=TRUE)
  update(CDF.model, burnin, progress.bar="none")
  CDF = jags.samples(CDF.model, variable.names=c("mu0","sigma0","pi0"), 
                     numb.iter,thin=numb.iter/5000, progress.bar="none")
  # CDF.coda = coda.samples(CDF.model, variable.names=c("mu0","sigma0","pi0"), 
  #                    numb.iter,thin=numb.iter/5000, progress.bar="none")
  
  #########################
  ##### JAGS estimates
  #########################
  mu0.hat <- apply(CDF$mu0,1:2,mean)
  sigma0.hat <- apply(CDF$sigma0,1:2,mean)
  pi0.hat <- apply(CDF$pi0,1:2,mean)
  G <- dim(mu0.hat)[2]
  
  ########################
  ###### Estimate PV
  ########################
  
  n1 = length(y1)
  eta.tab <-  rep(NA, n1)
  for(j in 1:n1)
  {
    eta.tab[j] <- 1-sum(apply(pi0.hat,1,mean)*pnorm(y1[j], apply(mu0.hat,1,mean), apply(sigma0.hat,1,mean)))
  }
  eta = eta.tab
  
  inverse = function(fn, interval = NULL, lower = min(interval), upper = max(interval), ...){
    Vectorize(function(y){
      uniroot(f=function(x){fn(x)-y}, lower=lower, upper=upper, ...)$root
    })
  }
  eta = ifelse(eta < 0.0000001, 0.0000001, 
               ifelse(eta > (1-0.0000001), (1-0.0000001), eta))
  
  #########################
  #### Stage 2
  #########################
  
  #########################
  #### FIt JAGS model
  #########################
  dat = list(n1=length(y1), H=30, z=qnorm(eta))
  PV.model = jags.model(textConnection(dp_PV), data = dat, n.chains = 3, quiet=TRUE)
  update(PV.model, burnin,progress.bar = "none")
  PVjags = jags.samples(PV.model, variable.names=c("mu","sigma","pi","a"), 
                        numb.iter,thin=numb.iter/5000,progress.bar = "none")
  
  #########################
  ##### JAGS estimates
  #########################
  mu.hat <- apply(PVjags$mu,1:2,mean)
  sigma.hat <- apply(PVjags$sigma,1:2,mean)
  pi.hat <- apply(PVjags$pi,1:2,mean)
  G <- dim(mu.hat)[2]

  ##################################
  ##### Estimate ROC 
  ##################################
  
  SPV.ROC.est = function(x){
    grid <- seq(0,1,length=100)
    ROC <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC[j] <- sum(pi.hat[,x]*pnorm(qnorm(grid[j]), mu.hat[,x], sigma.hat[,x]))
    }
    return(ROC)
  }
  SPV.AUC.est = function(x){
    grid <- seq(0,1,length=100)
    ROC <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC[j] <- sum(pi.hat[,x]*pnorm(qnorm(grid[j]), mu.hat[,x], sigma.hat[,x]))
    }
    AUC = trap(ROC)
    # AUC = ifelse(AUC < 0.5, 1-AUC, AUC)
    return(AUC)
  }
  
  ROC.tab = sapply(1:G, SPV.ROC.est)
  ROC.hat = apply(ROC.tab, 1, mean)
  AUC.tab = sapply(1:G, SPV.AUC.est)
  AUC.hat = c(mean(AUC.tab), sd(AUC.tab), quantile(AUC.tab, probs = c(0.025, 0.975)))
  AUC.hat = matrix(AUC.hat, nrow = 1, byrow = TRUE)
  colnames(AUC.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  SPV.est.func = function(x){
    
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
    
    grid <- seq(0,1,length=100)
    ROC <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC[j] <- sum(pi.hat[,x]*pnorm(qnorm(grid[j]), mu.hat[,x], sigma.hat[,x]))
    }
    
    AUC = trap(ROC)
    # AUC = ifelse(AUC < 0.5, 1-AUC, AUC)
    
    ##########################
    ###########
    ##### Optimal functions
    ###########
    ##########################
    
    sp = function(c){
      sum(pi0.hat[,x]*pnorm(c, mu0.hat[,x], sigma0.hat[,x]))
    }
    se = function(c){
      sum(pi.hat[,x]*pnorm(qnorm(1-sum(pi0.hat[,x]*pnorm(c, mu0.hat[,x], sigma0.hat[,x]))),
                           mu.hat[,x], sigma.hat[,x]))
    }
    

    J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
    ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
    CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
    IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
    
    J.est = c(J.c$par)
    ER.est = c(ER.c$par)
    CZ.est = c(CZ.c$par)
    IU.est = c(IU.c$par)
    cutoff.tab = c(J.est, ER.est, CZ.est, IU.est)
    # out = list(AUC = AUC, ROC = ROC, cutoff = cutoff.tab)
    # return(out)
    return(cutoff = cutoff.tab)
  }
  
  cutoff.all = sapply(1:G, SPV.est.func)
  # cutoff.all = -cutoff.all
  cutoff.hat = cbind(apply(cutoff.all, 1, mean),
                     apply(cutoff.all, 1, sd),
                     t(apply(cutoff.all, 1, function(x) quantile(x, probs = c(0.025, 0.975)))))
  rownames(cutoff.hat) = c("J","ER","CZ","IU")
  colnames(cutoff.hat) = c("Mean","SD","Low 95%", "Upper 95%")
  
  out=list(ROC=ROC.hat, AUC=AUC.hat, cutoff=cutoff.hat)
  return(out)
}


