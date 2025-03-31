

# rm(list=ls())
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

trap <- function(y, grid = seq(0,1,length=100)){
  # grid <- seq(0,1,length=100)
  val <- pracma::trapz(grid,y)
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

model.diagnostics.func = function(mcmc_samples){
  
  # Plot trace plots
  # plot(mcmc_samples)
  
  # Calculate Gelman-Rubin statistics
  gelman_rubin <- gelman.diag(mcmc_samples)
  # print(gelman_rubin)
  
  # Calculate autocorrelation
  autocorr <- autocorr.diag(mcmc_samples)
  # print(autocorr)
  
  # Calculate effective sample size
  eff_sample_size <- effectiveSize(mcmc_samples)
  # print(eff_sample_size)
  
  # Calculate Monte Carlo Standard Error (MCSE)
  n_eff <- eff_sample_size  # Effective sample size
  n_samples <- length(mcmc_samples[[1]])  # Total number of samples
  # Extract MCMC samples as a matrix
  mcmc_matrix <- as.matrix(mcmc_samples)
  
  # Calculate variance of the MCMC samples
  mcmc_var <- apply(mcmc_matrix, 2, function(x){var(x, na.rm = TRUE)})
  
  # Calculate MCSE
  mcse <- sqrt(mcmc_var / n_eff) * sqrt(n_samples / (n_samples - 1))
  # print(mcse)
  
  # Conduct convergence tests
  convergence_test <- gelman.plot(mcmc_samples)
  # print(convergence_test)
  
  out = list(gelman.rubin = gelman_rubin, autocorr = autocorr, eff_sample_size=eff_sample_size,
             mcse = mcse, convergence_test = convergence_test)
  return(out)
}


Binorm_diagnostics <- function(mu0,mu1,sigma0,sigma1, grid = seq(0,1,length=100)){
  b <- sigma0/sigma1
  a <- (mu1-mu0)/sigma1
  # grid <- seq(0,1,length=100)
  
  AUC <- pnorm(a/sqrt(1+b^2))
  ROC = pnorm(a + b*qnorm(grid))
  out=list(ROC = ROC, AUC = AUC,
           a = a, b = b)
  return(out)
}

Bigamma_diagnostics <- function(shape, scale0, scale1, grid = seq(0,1,length=100)){
  # grid <- seq(0,1,length=100)
  
  # shape = 0.5; scale0 = 0.1; scale1 = 0.2
  ROC = 1-pgamma(qgamma(1-grid, shape = shape, scale = scale0), shape = shape, scale = scale1)
  # AUC <- trap(ROC)
  AUC <- 1-pf(scale0/scale1, df1 = 2*shape, df2 = 2*shape)
  out=list(ROC = ROC, AUC = AUC)
  return(out)
}

Binorm_diagnostics_ab <- function(a,b, grid = seq(0,1,length=100)){
  # grid <- seq(0,1,length=100)
  AUC <- pnorm(a/sqrt(1+b^2))
  ROC = pnorm(a + b*qnorm(grid))
  out=list(ROC = ROC, AUC = AUC)
  return(out)
}


####################
####
####  Model fitting
####  function
####
####################

Empirical_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  # mn = mean(c(y0, y1))
  # s = sd(c(y0, y1))
  # y0 = (y0 - mn)/s
  # y1 = (y1 - mn)/s
  
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
  
  J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
  ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
  CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
  IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  # cutoff.tab[,1] = cutoff.tab[,1]*s + mn
  colnames(cutoff.tab) = c("cutoff","value","sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

NonPar_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  # mn = mean(c(y0, y1))
  # s = sd(c(y0, y1))
  # y0 = (y0 - mn)/s
  # y1 = (y1 - mn)/s
  
  n0 = length(y0)
  n1 = length(y1)
  h0 = 0.9*min(sd(y0), abs(diff(as.numeric(quantile(y0, probs = c(0.75, 0.25)))))/1.34 )*(n0^(-1.5))
  h1 = 0.9*min(sd(y1), abs(diff(as.numeric(quantile(y1, probs = c(0.75, 0.25)))))/1.34 )*(n1^(-1.5))
  
  F0 = function(c){
    CDF = (1/n0)*sum(pnorm((c-y0)/h0))
    return(CDF)
  }
  
  # inverse = function (f, lower = -100, upper = 100) {
  #   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  # }
  
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
  
  J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
  ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
  CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
  IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  # cutoff.tab[,1] = cutoff.tab[,1]*s + mn
  colnames(cutoff.tab) = c("cutoff","value","sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

BiNormal_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  # mn = mean(c(y0, y1))
  # s = sd(c(y0, y1))
  # y0 = (y0 - mn)/s
  # y1 = (y1 - mn)/s
  
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
  # binorm_1 = coda.samples(binorm.model, variable.names=c("mu0","mu1","sigma0","sigma1"), 
  # numb.iter,thin=numb.iter/5000, progress.bar="none")
  
  #MCMCsummary(binorm_1,params=c("mu0","mu1","sigma0","sigma1"),n.eff = TRUE)
  #MCMCplot(binorm_1,params=c("mu0","mu1","sigma0","sigma1"))
  #MCMCtrace(binorm_1,params=c("mu0","mu1","sigma0","sigma1"))
  
  #########################
  ##### JAGS estimates
  #########################
  mu0.hat <- apply(binorm$mu0,1:2,mean)
  mu1.hat <- apply(binorm$mu1,1:2,mean)
  sigma0.hat <- apply(binorm$sigma0,1:2,mean)
  sigma1.hat <- apply(binorm$sigma1,1:2,mean)
  
  mu0.hat = mean(mu0.hat)
  mu1.hat = mean(mu1.hat)
  sigma0.hat = mean(sigma0.hat)
  sigma1.hat = mean(sigma1.hat)
  
  a.hat = (mu1.hat-mu0.hat)/sigma1.hat
  b.hat = sigma0.hat/sigma1.hat
  grid = seq(0,1,length=100)
  
  # ROC <- array(NA,dim=c(length(grid),dim(a.hat)[2]))
  # for(i in 1:dim(a.hat)[2]){
  #   ROC[,i] <- pnorm(a.hat[i]+b.hat[i]*qnorm(grid))
  # }
  # AUC <- pnorm(a.hat/sqrt(1+b.hat^2))
  
  ROC <- pnorm(a.hat + b.hat*qnorm(grid))
  AUC <- pnorm(a.hat/sqrt(1+b.hat^2))
  
  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  se = function(c){
    1 - pnorm(c, mu1.hat, sigma1.hat)
    # 1 - pnorm(c,mean(mu1.hat), mean(sigma1.hat))
  }
  sp = function(c){
    pnorm(c, mu0.hat, sigma0.hat)
    # pnorm(c,mean(mu0.hat), mean(sigma0.hat))
  }
  
  J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
  ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
  CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
  IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
  
  J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
  ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
  CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
  IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  # cutoff.tab[,1] = cutoff.tab[,1]*s + mn
  colnames(cutoff.tab) = c("cutoff","value","sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

paramPV_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  # mn = mean(c(y0, y1))
  # s = sd(c(y0, y1))
  # y0 = (y0 - mn)/s
  # y1 = (y1 - mn)/s
  
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
  
  ##stage 1
  
  stage1.model = lm(y0~1)
  mu0 = as.numeric(coefficients(stage1.model))
  sigma0 = summary(stage1.model)$sigma
  
  PV <- 1-pnorm(y1,mu0, sigma0)
  PV <- ifelse(PV<(1e-15),(1e-15),ifelse(PV>(1-(1e-15)),(1-(1e-15)),PV))
  PV_one = PV
  
  #####################################
  #######    PV Model           #######
  #####################################
  
  stage2.model = lm(qnorm(PV_one)~1)
  mu = as.numeric(coefficients(stage2.model))
  sigma = summary(stage2.model)$sigma

  ROC <- pnorm(qnorm(grid),mu,sigma)
  AUC = trap(ROC)
  
  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  sp = function(c){
    # pnorm(c, mean(mu0), mean(sigma0))
    pnorm(c, mu0, sigma0)
  }
  se = function(c){
    # pnorm(1-pnorm(c, mean(mu0), mean(sigma0)), mean(mu), mean(sigma))
    # pnorm(1-pnorm(c, mu0, sigma0), mu, sigma)
    pnorm(qnorm(1-pnorm(c, mu0, sigma0)), mu, sigma)
  }
  
  # stage1.new = lm(y1~1)
  # mu1 = as.numeric(coefficients(stage1.new))
  # sigma1 = summary(stage1.new)$sigma
  # F1.regular = function(c){
  #   pnorm(c, mu1, sigma1)
  # }
  # F1.PV = function(c){
  #   1-pnorm(qnorm(1-pnorm(c,mu0,sigma0)), mu, sigma)
  # }
  
  J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
  ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
  CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
  IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
  
  J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
  ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
  CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
  IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  # cutoff.tab[,1] = cutoff.tab[,1]*s + mn
  colnames(cutoff.tab) = c("cutoff","value","sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

PV_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  # mn = mean(c(y0, y1))
  # s = sd(c(y0, y1))
  # y0 = (y0 - mn)/s
  # y1 = (y1 - mn)/s
  
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
  mu0 = mean(mu0)
  sigma0 = mean(sigma0)
  
  # PV <- array(NA, dim=c(length(y1),length(mu0)))
  # for(i in 1:length(mu0)){
  #   PV[,i] <- 1-pnorm(y1,mu0[i],sigma0[i])
  #   PV[,i] <- ifelse(PV[,i]<(1e-15),(1e-15),ifelse(PV[,i]>(1-(1e-15)),(1-(1e-15)),PV[,i]))
  # }
  # PV_one <- apply(PV,1,mean)
  
  PV <- 1-pnorm(y1,mu0, sigma0)
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
  mu = mean(mu)
  sigma = mean(sigma)

  # ROC <- array(NA, dim=c(length(grid),length(mu)))
  # AUC <- rep(NA,length(mu))
  # for (j in 1:length(mu)) {
  #   ROC[,j] <- pnorm(qnorm(grid),mu[j],sigma[j])
  #   AUC[j] <- trap(ROC[,j])
  # }
  
  ROC <- pnorm(qnorm(grid),mu,sigma)
  AUC = trap(ROC)

  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  sp = function(c){
    # pnorm(c, mean(mu0), mean(sigma0))
    pnorm(c, mu0, sigma0)
  }
  se = function(c){
    # pnorm(1-pnorm(c, mean(mu0), mean(sigma0)), mean(mu), mean(sigma))
    pnorm(qnorm(1-pnorm(c, mu0, sigma0)), mu, sigma)
  }
  
  J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
  ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
  CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
  IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
  
  J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
  ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
  CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
  IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  # cutoff.tab[,1] = cutoff.tab[,1]*s + mn
  colnames(cutoff.tab) = c("cutoff","value","sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
  return(out)
}

SemiPar_func <- function(y0, y1, burnin = 5000, numb.iter = 5000){
  
  # mn = mean(c(y0, y1))
  # s = sd(c(y0, y1))
  # y0 = (y0 - mn)/s
  # y1 = (y1 - mn)/s
  
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
    
    mu0.hat = apply(mu0.hat,1,mean)
    sigma0.hat = apply(sigma0.hat,1,mean)
    pi0.hat = apply(pi0.hat,1,mean)
    
    ########################
    ###### Estimate PV
    ########################
    
    # n1 = length(y1)
    # eta.tab <-  array(NA,dim=c(n1,G))
    # for(j in 1:n1)
    # {
    #   for(i in 1:G)
    #   {
    #     eta.tab[j,i] <- 1-sum(pi0.hat[,i]*pnorm((y1[j]-mu0.hat[,i])/sigma0.hat[,i]))
    #   }
    # }
    # eta=apply(eta.tab,1,mean)

    n1 = length(y1)
    eta.tab <-  rep(NA, n1)
    for(j in 1:n1)
    {
        # eta.tab[j] <- 1-sum(pi0.hat*pnorm((y1[j]-mu0.hat)/sigma0.hat))
        eta.tab[j] <- 1-sum(pi0.hat*pnorm(y1[j],mu0.hat,sigma0.hat))
    }
    eta = eta.tab
    
    inverse = function(fn, interval = NULL, lower = min(interval), upper = max(interval), ...){
      Vectorize(function(y){
        uniroot(f=function(x){fn(x)-y}, lower=lower, upper=upper, ...)$root
      })
    }
    F0 = function(c){
      sum(pi0.hat*pnorm(c, mu0.hat, sigma0.hat))
    }
    
    F0.inv = inverse(F0, lower=-100, upper=100)
    y1.new = F0.inv(1-eta)
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
    mu.hat = apply(mu.hat,1,mean)
    sigma.hat = apply(sigma.hat,1,mean)
    pi.hat = apply(pi.hat,1,mean)

    ##################################
    ##### Estimate CDF (ROC) at 
    ##### different iteration level
    ##################################
    grid <- seq(0,1,length=100)
    # ROC.tab <-  array(NA,dim=c(length(grid),G))
    # for(j in 1:length(grid))
    # {
    #   for(i in 1:G)
    #   {
    #     ROC.tab[j,i] <- sum(pi.hat[,i]*pnorm((logit(grid[j])-mu.hat[,i])/sigma.hat[,i]))
    #   }
    # }
    # ROC <- apply(ROC.tab,1,mean)
    
    ROC.tab <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
        # ROC.tab[j] <- sum(pi.hat*pnorm((logit(grid[j])-mu.hat)/sigma.hat))
      # ROC.tab[j] <- sum(pi.hat*pnorm(logit(grid[j]),mu.hat,sigma.hat))
      ROC.tab[j] <- sum(pi.hat*pnorm(qnorm(grid[j]),mu.hat,sigma.hat))
    }
    ROC = ROC.tab
    
    ##################
    ##### AUC
    ##################
    
    # AUC <- apply(ROC.tab,2,trap)
    AUC = trap(ROC)
    AUC = ifelse(AUC < 0.5, 1-AUC, AUC)
    
    ##########################
    ###########
    ##### Optimal functions
    ###########
    ##########################
    
    sp = function(c){
      # sum(apply(pi0.hat,1,mean)*pnorm(c, apply(mu0.hat,1,mean), apply(sigma0.hat,1,mean)))
      sum(pi0.hat*pnorm(c, mu0.hat, sigma0.hat))
    }
    se = function(c){
      # sum(apply(pi.hat,1,mean)*pnorm(1-sum(apply(pi0.hat,1,mean)*pnorm(c, apply(mu0.hat,1,mean), 
      #                                                              apply(sigma0.hat,1,mean))),
      #       apply(mu.hat,1,mean), apply(sigma.hat,1,mean)))
      sum(pi.hat*pnorm(qnorm(1-sum(pi0.hat*pnorm(c, mu0.hat, sigma0.hat))), mu.hat, sigma.hat))
    }
    
    # c.seq = seq(-3,3, length = 100)
    # plot(c.seq, pnorm(c.seq, 0.5, 0.9), type = "l")
    # lines(c.seq, sapply(c.seq, function(x){1-se(x)}), col="red")
    
    J.c = optim(mean(c(y0,y1)), J, method = "BFGS")
    ER.c = optim(mean(c(y0,y1)), ER, method = "BFGS")
    CZ.c = optim(mean(c(y0,y1)), CZ, method = "BFGS")
    IU.c = optim(mean(c(y0,y1)), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    # cutoff.tab[,1] = cutoff.tab[,1]*s + mn
    colnames(cutoff.tab) = c("cutoff","value","sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff.tab)
    return(out)
}

####################
####
####  Data generating
####  function
####
####################

BN.equal = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    mu0 = 0; mu1 = 0.2; sigma0 = 1; sigma1 = 1
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC
    ROC = bn$ROC
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, mu1, sigma1)
    
    se = function(c){
      1 - pnorm(c,mu1, sigma1)
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( ((mu0 + mu1)/2), J, method = "BFGS")
    ER.c = optim( ((mu0 + mu1)/2), ER, method = "BFGS")
    CZ.c = optim( ((mu0 + mu1)/2), CZ, method = "BFGS")
    IU.c = optim( ((mu0 + mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    mu0 = 0; mu1 = 1; sigma0 = 1; sigma1 = 1
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC
    ROC = bn$ROC
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, mu1, sigma1)
    
    se = function(c){
      1 - pnorm(c,mu1, sigma1)
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( ((mu0 + mu1)/2), J, method = "BFGS")
    ER.c = optim( ((mu0 + mu1)/2), ER, method = "BFGS")
    CZ.c = optim( ((mu0 + mu1)/2), CZ, method = "BFGS")
    IU.c = optim( ((mu0 + mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    mu0 = 0; mu1 = 2.5; sigma0 = 1; sigma1 = 1
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC
    ROC = bn$ROC
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, mu1, sigma1)
    
    se = function(c){
      1 - pnorm(c,mu1, sigma1)
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( ((mu0 + mu1)/2), J, method = "BFGS")
    ER.c = optim( ((mu0 + mu1)/2), ER, method = "BFGS")
    CZ.c = optim( ((mu0 + mu1)/2), CZ, method = "BFGS")
    IU.c = optim( ((mu0 + mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

BN.unequal = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    mu0 = 0; mu1 = 0.2; sigma0 = 1.2; sigma1 = 0.8
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, mu1, sigma1)
    
    se = function(c){
      1 - pnorm(c,mu1, sigma1)
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( ((mu0 + mu1)/2), J, method = "BFGS")
    ER.c = optim( ((mu0 + mu1)/2), ER, method = "BFGS")
    CZ.c = optim( ((mu0 + mu1)/2), CZ, method = "BFGS")
    IU.c = optim( ((mu0 + mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    mu0 = 0; mu1 = 1; sigma0 = 1.2; sigma1 = 0.5
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, mu1, sigma1)
    
    se = function(c){
      1 - pnorm(c,mu1, sigma1)
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( ((mu0 + mu1)/2), J, method = "BFGS")
    ER.c = optim( ((mu0 + mu1)/2), ER, method = "BFGS")
    CZ.c = optim( ((mu0 + mu1)/2), CZ, method = "BFGS")
    IU.c = optim( ((mu0 + mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    mu0 = 1; mu1 = 2.9; sigma0 = 0.5; sigma1 = 1.2
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, mu1, sigma1)
    
    se = function(c){
      1 - pnorm(c,mu1, sigma1)
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( ((mu0 + mu1)/2), J, method = "BFGS")
    ER.c = optim( ((mu0 + mu1)/2), ER, method = "BFGS")
    CZ.c = optim( ((mu0 + mu1)/2), CZ, method = "BFGS")
    IU.c = optim( ((mu0 + mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

Skewed.I = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    # mu0 = 2.5; mu1 = 2.55; sigma0 = 0.09; sigma1 = 0.25
    mu0 = 0; mu1 = 0.2; sigma0 = 1; sigma1 = 1; 
    
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = (rnorm(N, mu0, sigma0))^(2)
    y1 = (rnorm(N, mu1, sigma1))^(2)
    
    se = function(c){
      1 - pnorm(c^(1/2),mu1, sigma1)
      # 1 - ecdf(y1)(c)
    }
    sp = function(c){
      pnorm(c^(1/2),mu0, sigma0)
      # ecdf(y0)(c)
    }
    
    J.c = optim((mu0+mu1)^(2), J, method = "BFGS")
    ER.c = optim((mu0+mu1)^(2), ER, method = "BFGS")
    CZ.c = optim((mu0+mu1)^(2), CZ, method = "BFGS")
    IU.c = optim((mu0+mu1)^(2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    mu0 = 0; mu1 = 1; sigma0 = 1; sigma1 = 0.7;
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = (rnorm(N, mu0, sigma0))^(2)
    y1 = (rnorm(N, mu1, sigma1))^(2)
    
    se = function(c){
      1 - pnorm(c^(1/2),mu1, sigma1)
      # 1 - ecdf(y1)(c)
    }
    sp = function(c){
      pnorm(c^(1/2),mu0, sigma0)
      # ecdf(y0)(c)
    }
    
    J.c = optim((mu0+mu1)^(2), J, method = "BFGS")
    ER.c = optim((mu0+mu1)^(2), ER, method = "BFGS")
    CZ.c = optim((mu0+mu1)^(2), CZ, method = "BFGS")
    IU.c = optim((mu0+mu1)^(2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    mu0 = 1; mu1 = 2.5; sigma0 = 1; sigma1 = 0.5;
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = (rnorm(N, mu0, sigma0))^(2)
    y1 = (rnorm(N, mu1, sigma1))^(2)
    
    se = function(c){
      1 - pnorm(c^(1/2),mu1, sigma1)
      # 1 - ecdf(y1)(c)
    }
    sp = function(c){
      pnorm(c^(1/2),mu0, sigma0)
      # ecdf(y0)(c)
    }
    
    J.c = optim((mu0+mu1)^(2), J, method = "BFGS")
    ER.c = optim((mu0+mu1)^(2), ER, method = "BFGS")
    CZ.c = optim((mu0+mu1)^(2), CZ, method = "BFGS")
    IU.c = optim((mu0+mu1)^(2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

Skewed.II = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    # mu0 = 2.5; mu1 = 2.55; sigma0 = 0.09; sigma1 = 0.25
    mu0 = 0; mu1 = 0.2; sigma0 = 1; sigma1 = 1;
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = exp(rnorm(N, mu0, sigma0))
    y1 = exp(rnorm(N, mu1, sigma1))
    
    se = function(c){
      1 - pnorm(log(c),mu1, sigma1)
      #  1 - ecdf(y1)(c)
    }
    sp = function(c){
      pnorm(log(c),mu0, sigma0)
      # ecdf(y0)(c)
    }
    
    J.c = optim(exp((mu0+mu1)/2), J, method = "BFGS")
    ER.c = optim(exp((mu0+mu1)/2), ER, method = "BFGS")
    CZ.c = optim(exp((mu0+mu1)/2), CZ, method = "BFGS")
    IU.c = optim(exp((mu0+mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    mu0 = 0; mu1 = 1; sigma0 = 1; sigma1 = 0.7;
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = exp(rnorm(N, mu0, sigma0))
    y1 = exp(rnorm(N, mu1, sigma1))
    
    se = function(c){
      1 - pnorm(log(c),mu1, sigma1)
      #  1 - ecdf(y1)(c)
    }
    sp = function(c){
      pnorm(log(c),mu0, sigma0)
      # ecdf(y0)(c)
    }
    
    J.c = optim(exp((mu0+mu1)/2), J, method = "BFGS")
    ER.c = optim(exp((mu0+mu1)/2), ER, method = "BFGS")
    CZ.c = optim(exp((mu0+mu1)/2), CZ, method = "BFGS")
    IU.c = optim(exp((mu0+mu1)/2), IU, method = "BFGS")

    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    mu0 = 1; mu1 = 2.5; sigma0 = 1; sigma1 = 0.5;
    grid = seq(0,1,length=100)
    bn = Binorm_diagnostics(mu0,mu1,sigma0,sigma1,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = exp(rnorm(N, mu0, sigma0))
    y1 = exp(rnorm(N, mu1, sigma1))
    
    se = function(c){
      1 - pnorm(log(c),mu1, sigma1)
      #  1 - ecdf(y1)(c)
    }
    sp = function(c){
      pnorm(log(c),mu0, sigma0)
      # ecdf(y0)(c)
    }
    
    J.c = optim(exp((mu0+mu1)/2), J, method = "BFGS")
    ER.c = optim(exp((mu0+mu1)/2), ER, method = "BFGS")
    CZ.c = optim(exp((mu0+mu1)/2), CZ, method = "BFGS")
    IU.c = optim(exp((mu0+mu1)/2), IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

Skewed.III = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    shape = 0.5; scale0 = 0.1; scale1 = 0.15
    grid = seq(0,1,length=100)
    bg = Bigamma_diagnostics(shape, scale0, scale1, grid)
    AUC = bg$AUC; #AUC
    ROC = bg$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rgamma(N, shape = shape, scale = scale0)
    y1 = rgamma(N, shape = shape, scale = scale1)
    
    se = function(c){
       1 - pgamma(c, shape = shape, scale = scale1)
    }
    sp = function(c){
      pgamma(c, shape = shape, scale = scale0)
    }
    
    J.c = optim( (shape*(scale0 + scale1))/2, J, method = "BFGS")
    ER.c = optim( (shape*(scale0 + scale1))/2, ER, method = "BFGS")
    CZ.c = optim( (shape*(scale0 + scale1))/2, CZ, method = "BFGS")
    IU.c = optim( (shape*(scale0 + scale1))/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    shape = 0.5; scale0 = 0.1; scale1 = 0.6
    grid = seq(0,1,length=100)
    bg = Bigamma_diagnostics(shape, scale0, scale1, grid)
    AUC = bg$AUC; #AUC
    ROC = bg$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rgamma(N, shape = shape, scale = scale0)
    y1 = rgamma(N, shape = shape, scale = scale1)
    
    se = function(c){
      1 - pgamma(c, shape = shape, scale = scale1)
    }
    sp = function(c){
      pgamma(c, shape = shape, scale = scale0)
    }
    
    J.c = optim( (shape*(scale0 + scale1))/2, J, method = "BFGS")
    ER.c = optim( (shape*(scale0 + scale1))/2, ER, method = "BFGS")
    CZ.c = optim( (shape*(scale0 + scale1))/2, CZ, method = "BFGS")
    IU.c = optim( (shape*(scale0 + scale1))/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    shape = 0.5; scale0 = 0.1; scale1 = 7
    grid = seq(0,1,length=100)
    bg = Bigamma_diagnostics(shape, scale0, scale1, grid)
    AUC = bg$AUC; #AUC
    ROC = bg$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rgamma(N, shape = shape, scale = scale0)
    y1 = rgamma(N, shape = shape, scale = scale1)
    
    se = function(c){
      1 - pgamma(c, shape = shape, scale = scale1)
    }
    sp = function(c){
      pgamma(c, shape = shape, scale = scale0)
    }
    
    J.c = optim( (shape*(scale0 + scale1))/2, J, method = "BFGS")
    ER.c = optim( (shape*(scale0 + scale1))/2, ER, method = "BFGS")
    CZ.c = optim( (shape*(scale0 + scale1))/2, CZ, method = "BFGS")
    IU.c = optim( (shape*(scale0 + scale1))/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

Mixed.I = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    mu0 = 0; sigma0 = 1;
    mu11 = 0; sigma11 = 1;
    mu12 = 1; sigma12 = 5;
    pi = 0.5;
    grid = seq(0,1,length=100)
    
    a = (pi*mu11 + (1-pi)*mu12)/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    b = sigma0/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    
    bn = Binorm_diagnostics_ab(a,b,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, pi*mu11 + (1-pi)*mu12, sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2))
    
    se = function(c){
      1 - pnorm(c, pi*mu11 + (1-pi)*mu12, sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2))
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, J, method = "BFGS")
    ER.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, ER, method = "BFGS")
    CZ.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, CZ, method = "BFGS")
    IU.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    mu0 = 0; sigma0 = 1;
    mu11 = 0; sigma11 = 1;
    mu12 = 4; sigma12 = 5;
    pi = 0.5;
    grid = seq(0,1,length=100)
    
    a = (pi*mu11 + (1-pi)*mu12)/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    b = sigma0/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    
    bn = Binorm_diagnostics_ab(a,b,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, pi*mu11 + (1-pi)*mu12, sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2))
    
    se = function(c){
      1 - pnorm(c, pi*mu11 + (1-pi)*mu12, sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2))
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, J, method = "BFGS")
    ER.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, ER, method = "BFGS")
    CZ.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, CZ, method = "BFGS")
    IU.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    mu0 = 0; sigma0 = 1;
    mu11 = 0; sigma11 = 1;
    mu12 = 8; sigma12 = 5;
    pi = 0.5;
    grid = seq(0,1,length=100)
    
    a = (pi*mu11 + (1-pi)*mu12)/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    b = sigma0/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    
    bn = Binorm_diagnostics_ab(a,b,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, mu0, sigma0)
    y1 = rnorm(N, pi*mu11 + (1-pi)*mu12, sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2))
    
    se = function(c){
      1 - pnorm(c, pi*mu11 + (1-pi)*mu12, sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2))
    }
    sp = function(c){
      pnorm(c,mu0, sigma0)
    }
    
    J.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, J, method = "BFGS")
    ER.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, ER, method = "BFGS")
    CZ.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, CZ, method = "BFGS")
    IU.c = optim( (pi*mu11 + (1-pi)*mu12 + mu0)/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

Mixed.II = function(seed, N=50, sep = "Low"){
  
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
  
  if(sep == "Low"){
    mu01 = 0; sigma01 = 1;
    mu02 = 1; sigma02 = 2;
    pi0 = 0.5;
    mu11 = 0; sigma11 = 1;
    mu12 = 1.5; sigma12 = 2.5;
    pi1 = 0.4;
    grid = seq(0,1,length=100)
    
    a = ((pi1*mu11 + (1-pi1)*mu12)-(pi0*mu01 + (1-pi0)*mu02))/sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2)
    b = sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2)/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    
    bn = Binorm_diagnostics_ab(a,b,grid)
    AUC = bn$AUC; AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, pi0*mu01 + (1-pi0)*mu02, sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2))
    y1 = rnorm(N, pi1*mu11 + (1-pi1)*mu12, sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2))
    
    # plot(density(y0))
    # lines(density(y1))
    
    se = function(c){
      1 - pnorm(c, pi1*mu11 + (1-pi1)*mu12, sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2))
    }
    sp = function(c){
      pnorm(c, pi0*mu01 + (1-pi0)*mu02, sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2))
    }
    
    J.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, J, method = "BFGS")
    ER.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, ER, method = "BFGS")
    CZ.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, CZ, method = "BFGS")
    IU.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "Medium"){
    mu01 = 0; sigma01 = 1;
    mu02 = 1; sigma02 = 2;
    pi0 = 0.5;
    mu11 = 0; sigma11 = 1;
    mu12 = 2.5; sigma12 = 2.5;
    pi1 = 0.4;
    grid = seq(0,1,length=100)
    
    a = ((pi1*mu11 + (1-pi1)*mu12)-(pi0*mu01 + (1-pi0)*mu02))/sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2)
    b = sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2)/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    
    bn = Binorm_diagnostics_ab(a,b,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, pi0*mu01 + (1-pi0)*mu02, sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2))
    y1 = rnorm(N, pi1*mu11 + (1-pi1)*mu12, sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2))
    
    se = function(c){
      1 - pnorm(c, pi1*mu11 + (1-pi1)*mu12, sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2))
    }
    sp = function(c){
      pnorm(c, pi0*mu01 + (1-pi0)*mu02, sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2))
    }
    
    J.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, J, method = "BFGS")
    ER.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, ER, method = "BFGS")
    CZ.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, CZ, method = "BFGS")
    IU.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  else if(sep == "High"){
    mu01 = 0; sigma01 = 1;
    mu02 = 1; sigma02 = 2;
    pi0 = 0.5;
    mu11 = 0; sigma11 = 1;
    mu12 = 5; sigma12 = 2.5;
    pi1 = 0.4;
    grid = seq(0,1,length=100)
    
    a = ((pi1*mu11 + (1-pi1)*mu12)-(pi0*mu01 + (1-pi0)*mu02))/sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2)
    b = sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2)/sqrt((pi^2)*sigma11^2 + ((1-pi)^2)*sigma12^2)
    
    bn = Binorm_diagnostics_ab(a,b,grid)
    AUC = bn$AUC; #AUC
    ROC = bn$ROC; #plot(grid, ROC, type="l")
    
    set.seed(seed)
    y0 = rnorm(N, pi0*mu01 + (1-pi0)*mu02, sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2))
    y1 = rnorm(N, pi1*mu11 + (1-pi1)*mu12, sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2))
    
    se = function(c){
      1 - pnorm(c, pi1*mu11 + (1-pi1)*mu12, sqrt((pi1^2)*sigma11^2 + ((1-pi1)^2)*sigma12^2))
    }
    sp = function(c){
      pnorm(c, pi0*mu01 + (1-pi0)*mu02, sqrt((pi0^2)*sigma01^2 + ((1-pi0)^2)*sigma02^2))
    }
    
    J.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, J, method = "BFGS")
    ER.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, ER, method = "BFGS")
    CZ.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, CZ, method = "BFGS")
    IU.c = optim(((pi1*mu11 + (1-pi1)*mu12)+(pi0*mu01 + (1-pi0)*mu02))/2, IU, method = "BFGS")
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par), sp(J.c$par), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par), sp(ER.c$par), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par), sp(CZ.c$par), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par), sp(IU.c$par), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    # rbind(J.est, ER.est, CZ.est, IU.est)
  }
  
  cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
  colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                           "sensitivity (data)", "specificity (data)")
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff.tab)
  return(out)
}

####################
####
####  All fit
####
####################

cutoff.func = function(seed, N = 100, data.gen = BN.equal, sep = "Medium",
                       burnin = 5000, numb.iter = 5000){
  
  # seed = 2000; N = 100; data.gen = Mixed.II; sep = "Low";
  # burnin = 5000; numb.iter = 5000
  
  dat.mod = do.call(data.gen, list(seed, N, sep))
  
  y0 = dat.mod$y0; y1 = dat.mod$y1
  grid = dat.mod$grid
  TrueAUC = dat.mod$TrueAUC
  TrueROC = dat.mod$TrueROC
  True.cutoff = dat.mod$cutoff.tab
  
  mn = mean(c(y0, y1))
  s = sd(c(y0, y1))
  y0 = (y0 - mn)/s
  y1 = (y1 - mn)/s
  
  # tt = Sys.time()
  Emp.fit = Empirical_func(y0, y1, burnin = burnin, numb.iter = numb.iter)
  BN.fit = BiNormal_func(y0, y1, burnin = burnin, numb.iter = numb.iter)
  NonPar.fit = NonPar_func(y0, y1, burnin = burnin, numb.iter = numb.iter)
  PV.fit = PV_func(y0, y1, burnin = burnin, numb.iter = numb.iter)
  SemiPV.fit = SemiPar_func(y0, y1, burnin = burnin, numb.iter = numb.iter)
  # Sys.time() - tt
  
  AUC = data.frame(Emp = Emp.fit$AUC, BN = BN.fit$AUC, 
                   NonPar = NonPar.fit$AUC,
             PV = PV.fit$AUC, Semi.PV = SemiPV.fit$AUC)
  ROC = data.frame(grid = grid, Emp = Emp.fit$ROC, BN = BN.fit$ROC, 
                   NonPar = NonPar.fit$ROC,
                   PV = PV.fit$ROC, Semi.PV = SemiPV.fit$ROC)
  
  Emp.fit$cutoff.tab[,1] = Emp.fit$cutoff.tab[,1]*s + mn
  BN.fit$cutoff.tab[,1] = BN.fit$cutoff.tab[,1]*s + mn
  NonPar.fit$cutoff.tab[,1] = NonPar.fit$cutoff.tab[,1]*s + mn
  PV.fit$cutoff.tab[,1] = PV.fit$cutoff.tab[,1]*s + mn
  SemiPV.fit$cutoff.tab[,1] = SemiPV.fit$cutoff.tab[,1]*s + mn

  J.tab = data.frame(Emp = Emp.fit$cutoff.tab[1,], BN = BN.fit$cutoff.tab[1,], 
                     NonPar = NonPar.fit$cutoff.tab[1,],
                     PV = PV.fit$cutoff.tab[1,], 
                     Semi.PV = SemiPV.fit$cutoff.tab[1,])
  ER.tab = data.frame(Emp = Emp.fit$cutoff.tab[2,], BN = BN.fit$cutoff.tab[2,], 
                      NonPar = NonPar.fit$cutoff.tab[2,],
                      PV = PV.fit$cutoff.tab[2,], 
                      Semi.PV = SemiPV.fit$cutoff.tab[2,])
  CZ.tab = data.frame(Emp = Emp.fit$cutoff.tab[3,], BN = BN.fit$cutoff.tab[3,], 
                      NonPar = NonPar.fit$cutoff.tab[3,],
                      PV = PV.fit$cutoff.tab[3,], 
                      Semi.PV = SemiPV.fit$cutoff.tab[3,])
  IU.tab = data.frame(Emp = Emp.fit$cutoff.tab[4,], BN = BN.fit$cutoff.tab[4,], 
                      NonPar = NonPar.fit$cutoff.tab[4,],
                      PV = PV.fit$cutoff.tab[4,], 
                      Semi.PV = SemiPV.fit$cutoff.tab[4,])
  
  out = list(grid = grid, TrueAUC = TrueAUC, TrueROC = TrueROC, True.cutoff = True.cutoff,
             AUC = AUC, ROC = ROC, J.tab = J.tab, ER.tab = ER.tab,
             CZ.tab = CZ.tab, IU.tab = IU.tab)
  return(out)
}


final.fit = function(seed, n, sep){
  
  # source("Nocov sim func.R")
  
  ## possible variations of sep: "Low", "Medium", "High"
  ## possible variations of n: 50, 100, 500
  ## possible variations of seed: any number
  
  fit.BN.equal = cutoff.func(seed = seed, N = n, data.gen = BN.equal, sep = sep)
  fit.BN.unequal = cutoff.func(seed = seed, N = n, data.gen = BN.unequal, sep = sep)
  fit.Skewed.I = cutoff.func(seed = seed, N = n, data.gen = Skewed.I, sep = sep)
  fit.Skewed.II = cutoff.func(seed = seed, N = n, data.gen = Skewed.II, sep = sep)
  fit.Skewed.III = cutoff.func(seed = seed, N = n, data.gen = Skewed.III, sep = sep)
  fit.Mixed.I = cutoff.func(seed = seed, N = n, data.gen = Mixed.I, sep = sep)
  fit.Mixed.II = cutoff.func(seed = seed, N = n, data.gen = Mixed.II, sep = sep)
  
  out = list(fit.BN.equal = fit.BN.equal, fit.BN.unequal = fit.BN.unequal, fit.Skewed.I = fit.Skewed.I,
             fit.Skewed.II = fit.Skewed.II, fit.Skewed.III = fit.Skewed.III, fit.Mixed.I = fit.Mixed.I,
             fit.Mixed.II = fit.Mixed.II)
  
  return(out)
}



