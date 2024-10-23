
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

Binorm_diagnostics_ab <- function(a,b, grid = seq(0,1,length=100)){
  # grid <- seq(0,1,length=100)
  AUC <- pnorm(a/sqrt(1+b^2))
  ROC = pnorm(a + b*qnorm(grid))
  out=list(ROC = ROC, AUC = AUC)
  return(out)
}

Binorm_diagnostics_linear_cov <- function(a01,a11,a02,a12,sigma0,sigma1,cov, 
                                          grid = seq(0,1,length=100)){
  b <- sigma0/sigma1
  a <- ((a02-a01)+(a12-a11)*cov)/sigma1
  # grid = seq(0,1,length=100)
  
  AUC <- pnorm(a/sqrt(1+b^2))
  ROC = pnorm(a + b*qnorm(grid))
  out=list(ROC = ROC, AUC = AUC,
           a = a, b = b)
  return(out)
}

Binorm_diagnostics_mixed <- function(a00, a01, a101, a111, a102, a112, pi,
                                     sigma0,sigma1,cov, 
                                     grid = seq(0,1,length=100)){
  
  b <- sigma0/sigma1
  a = (pi*(a101-a102)+(a102-a00)+(a112-a01)*cov+pi*(a111-a112)*cov)/sigma1
  # grid = seq(0,1,length=100)
  
  AUC <- pnorm(a/sqrt(1+b^2))
  ROC = pnorm(a + b*qnorm(grid))
  out=list(ROC = ROC, AUC = AUC,
           a = a, b = b)
  return(out)
}


####################
####
####  Data generating
####  function
####
####################

BN = function(seed, N=50, X0, X1, c = c(0,1)){
  
  b00 = 1; b01 = 1
  b10 = 1.5; b11 = 2
  sd0 = sd1 = 1
  bn0 = Binorm_diagnostics_linear_cov(a01=b00,a11=b01,
                                     a02=b10,a12=b11,
                                     sigma0=sd0,sigma1=sd1,cov = c[1])
  bn1 = Binorm_diagnostics_linear_cov(a01=b00,a11=b01,
                                      a02=b10,a12=b11,
                                      sigma0=sd0,sigma1=sd1,cov = c[2])
  AUC = c(bn0$AUC, bn1$AUC #, bn2$AUC
          )
  ROC = cbind(bn0$ROC, bn1$ROC #, bn2$ROC
              )
  grid = seq(0,1,length = 100)
  
  set.seed(seed)
  y0 = rnorm(N,b00+b01*X0, sd0)
  y1 = rnorm(N,b10+b11*X1, sd1)

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
  
  se = function(c,x){
    1 - pnorm(c,b10+b11*x, sd1)
  }
  sp = function(c,x){
    pnorm(c,b00+b01*x, sd0)
  }
  
  summary.cutoff = function(c){
    J.c = optim(par = mean(c(b00+b01*c,b10+b11*c)), fn = J, method = "BFGS", x = c)
    ER.c = optim(par = mean(c(b00+b01*c,b10+b11*c)), fn = ER, method = "BFGS", x = c)
    CZ.c = optim(par = mean(c(b00+b01*c,b10+b11*c)), fn = CZ, method = "BFGS", x = c)
    IU.c = optim(par = mean(c(b00+b01*c,b10+b11*c)), fn = IU, method = "BFGS", x = c)
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par, x = c), sp(J.c$par, x = c), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = c), sp(ER.c$par, x = c), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = c), sp(CZ.c$par, x = c), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = c), sp(IU.c$par, x = c), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    return(cutoff.tab)
  }
  
  cutoff = list(
  cutoff0 = summary.cutoff(c = c[1]),
  cutoff1 = summary.cutoff(c = c[2]) #,
  # cutoff2 = summary.cutoff(c = c[3])
  )
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff)
  return(out)
}
Skewed = function(seed, N=50, X0, X1, c = c(0,1)){
  
  a00 = 3; a01 = 0.1
  a10 = 5; a11 = 9
  shape0 = 2
  shape1 = 2

  grid = seq(0,1,length=100)
  bg0 = Bigamma_diagnostics_linear_cov(a00=a00,a01=a01,
                                       a10=a10,a11=a11,
                                      shape0=shape0,shape1=shape1,cov = c[1])
  bg1 = Bigamma_diagnostics_linear_cov(a00=a00,a01=a01,
                                       a10=a10,a11=a11,
                                      shape0=shape0,shape1=shape1,cov = c[2])

  AUC = c(bg0$AUC, bg1$AUC #, bg2$AUC
          ); AUC
  ROC = cbind(bg0$ROC, bg1$ROC #, bg2$ROC
              )
  grid = seq(0,1,length = 100)
  
  set.seed(seed)
  y0 = rgamma(N, shape = shape0, scale = (a00+a01*X0))
  y1 = rgamma(N, shape = shape1, scale = (a10+a11*X1))
  
  
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
  
  se = function(c,x){
    # 1 - pgamma(c, shape = shape1, scale = exp(b10+b11*x))
    1 - pgamma(c, shape = shape1, scale = (a10+a11*x))
  }
  sp = function(c,x){
    # pgamma(c,shape = shape0, scale = exp(b00+b01*x))
    pgamma(c,shape = shape0, scale = (a00+a01*x))
  }
  
  summary.cutoff = function(c){
    J.c = optim(par = mean(c(shape0*(a00+a01*c),shape1*(a10+a11*c))), fn = J, method = "BFGS", x = c)
    ER.c = optim(par = mean(c(shape0*(a00+a01*c),shape1*(a10+a11*c))), fn = ER, method = "BFGS", x = c)
    CZ.c = optim(par = mean(c(shape0*(a00+a01*c),shape1*(a10+a11*c))), fn = CZ, method = "BFGS", x = c)
    IU.c = optim(par = mean(c(shape0*(a00+a01*c),shape1*(a10+a11*c))), fn = IU, method = "BFGS", x = c)
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par, x = c), sp(J.c$par, x = c), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = c), sp(ER.c$par, x = c), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = c), sp(CZ.c$par, x = c), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = c), sp(IU.c$par, x = c), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    return(cutoff.tab)
  }
  
  cutoff = list(
    cutoff0 = summary.cutoff(c = c[1]),
    cutoff1 = summary.cutoff(c = c[2]) #,
    # cutoff2 = summary.cutoff(c = c[3])
  )
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff)
  return(out)
}
Mixed = function(seed, N=50, X0, X1, c = c(0,1)){
  
  a00 = 0; a01 = 1
  a101 = 0; a111 = 1
  a102 = 1; a112 = 5
  pi = 0.5
  # X0 = runif(N)
  # X1 = runif(N,1,3)
  sd0 = 1
  sd1 = 1.5
  
  grid = seq(0,1,length=100)
  
  bm0 = Binorm_diagnostics_mixed(a00, a01, a101, a111, a102, a112, pi,
                                       sigma0 = sd0,sigma1=sd1,cov = c[1])
  bm1 = Binorm_diagnostics_mixed(a00, a01, a101, a111, a102, a112, pi,
                                 sigma0 = sd0,sigma1=sd1,cov = c[2])

  AUC = c(bm0$AUC, bm1$AUC #, bm2$AUC
          ); AUC
  ROC = cbind(bm0$ROC, bm1$ROC #, bm2$ROC
              )
  grid = seq(0,1,length = 100)
  
  set.seed(seed)
  mu0 = a00 + a01*X0
  mu11 = a101 + a111*X1
  mu12 = a102 + a112*X1
  y0 = rnorm(N, mu0, sd0)
  y1 = rnorm(N, pi*mu11 + (1-pi)*mu12, sd1)
  
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
  
  se = function(c,x){
    1 - pnorm(c, mean = pi*(a101+a111*x) + (1-pi)*(a102+a112*x), sd = sd1)
  }
  sp = function(c,x){
    pnorm(c, mean = a00+a01*x, sd = sd0)
  }
  
  summary.cutoff = function(c){
    J.c = optim(par = mean(c(a00 + a01*c, pi*(a101 + a111*c) + (1-pi)*(a102 + a112*c))), fn = J, method = "BFGS", x = c)
    ER.c = optim(par = mean(c(a00 + a01*c, pi*(a101 + a111*c) + (1-pi)*(a102 + a112*c))), fn = ER, method = "BFGS", x = c)
    CZ.c = optim(par = mean(c(a00 + a01*c, pi*(a101 + a111*c) + (1-pi)*(a102 + a112*c))), fn = CZ, method = "BFGS", x = c)
    IU.c = optim(par = mean(c(a00 + a01*c, pi*(a101 + a111*c) + (1-pi)*(a102 + a112*c))), fn = IU, method = "BFGS", x = c)
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par, x = c), sp(J.c$par, x = c), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = c), sp(ER.c$par, x = c), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = c), sp(CZ.c$par, x = c), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = c), sp(IU.c$par, x = c), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    return(cutoff.tab)
  }
  
  cutoff = list(
    cutoff0 = summary.cutoff(c = c[1]),
    cutoff1 = summary.cutoff(c = c[2]) #,
    # cutoff2 = summary.cutoff(c = c[3])
  )
  
  out = list(y0 = y0, y1 = y1, TrueAUC = AUC, 
             grid = grid, TrueROC = ROC,
             cutoff.tab = cutoff)
  return(out)
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
  
  b00.hat = mean(b00.hat)
  b10.hat = mean(b10.hat)
  b01.hat = mean(b01.hat)
  b11.hat = mean(b11.hat)
  sigma0.hat = mean(sigma0.hat)
  sigma1.hat = mean(sigma1.hat)
  
  est.tab = data.frame(b00 = b00.hat, b10 = b10.hat,
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
  
  bn0 = BN_accuracy_cutoff(y0, y1, b00=b00, b01=b01, b10=b10, b11=b11, 
                     sigma0=sigma0, sigma1=sigma1, x=x[1])
  bn1 = BN_accuracy_cutoff(y0, y1, b00=b00, b01=b01, b10=b10, b11=b11, 
                           sigma0=sigma0, sigma1=sigma1, x=x[2])

  AUC = c(bn0$AUC, bn1$AUC #, bn3$AUC
          )
  ROC = cbind(bn0$ROC, bn1$ROC #, bn3$ROC
              )
  cutoff = list(cutoff0 = bn0$cutoff.tab,
                cutoff1 = bn1$cutoff.tab #,
                # cutoff2 = bn3$cutoff.tab
                )
  out = list(AUC = AUC, ROC = ROC, cutoff = cutoff)
  return(out)
}


freqPV_reg_fit <- function(y0, y1, X0, X1, x){
  
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
  
  grid <- seq(0,1,length=100)
  
  ##stage 1
  
  stage1.model = lm(y0~X0)
  b00.hat = as.numeric(coefficients(stage1.model))[1]
  b01.hat = as.numeric(coefficients(stage1.model))[2]
  sigma0 = summary(stage1.model)$sigma
  
  PV <- 1-pnorm(y1, b00.hat + b01.hat*X1, sigma0)
  PV <- ifelse(PV<(1e-15),(1e-15),ifelse(PV>(1-(1e-15)),(1-(1e-15)),PV))
  PV_one = PV
  
  #####################################
  #######    PV Model           #######
  #####################################
  
  stage2.model = lm(qnorm(PV_one)~X1)
  b0 = as.numeric(coefficients(stage2.model))[1]
  b1 = as.numeric(coefficients(stage2.model))[2]
  sigma = summary(stage2.model)$sigma
  
  ROC = sapply(x, function(c){pnorm(qnorm(grid), b0+b1*c, sigma)})
  AUC = apply(ROC,2,trap)

  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  se = function(c, x){
    pnorm(qnorm(1 - pnorm(c, b00.hat + b01.hat*x, sigma0)), b0 + b1*x, sigma)
  }
  
  sp = function(c, x){
    pnorm(c, b00.hat + b01.hat*x, sigma0)
  }

  summary.cutoff = function(c){
    J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS", x = c)
    ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS", x = c)
    CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS", x = c)
    IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS", x = c)
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par, x = c), sp(J.c$par, x = c), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = c), sp(ER.c$par, x = c), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = c), sp(CZ.c$par, x = c), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = c), sp(IU.c$par, x = c), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    return(cutoff.tab)
  }
  
  cutoff = list(
    cutoff0 = summary.cutoff(c = x[1]),
    cutoff1 = summary.cutoff(c = x[2]) #,
    # cutoff3 = summary.cutoff(c = x[3])
  )
  
  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff)
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
  sigma0 = mean(sigma0)
  
  PV <- 1-pnorm(y1, b00.hat + b01.hat*X1, sigma0)
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
  b0.hat = mean(b0)
  b1.hat = mean(b1)
  sigma.hat = mean(sigma)
  
  ROC = sapply(x, function(c){pnorm(qnorm(grid), b0.hat + b1.hat*c, sigma.hat)})
  AUC = apply(ROC,2,trap)
  
  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  se = function(c, x){
    pnorm(qnorm(1 - pnorm(c, b00.hat + b01.hat*x, sigma0)), b0.hat + b1.hat*x, sigma.hat)
  }
  
  sp = function(c, x){
    pnorm(c, b00.hat + b01.hat*x, sigma0)
  }
  
  summary.cutoff = function(c){
    J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS", x = c)
    ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS", x = c)
    CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS", x = c)
    IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS", x = c)
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par, x = c), sp(J.c$par, x = c), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = c), sp(ER.c$par, x = c), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = c), sp(CZ.c$par, x = c), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = c), sp(IU.c$par, x = c), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    return(cutoff.tab)
  }
  
  cutoff = list(
    cutoff0 = summary.cutoff(c = x[1]),
    cutoff1 = summary.cutoff(c = x[2]) #,
    # cutoff3 = summary.cutoff(c = x[3])
  )

  out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff)
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
  # CDF.coda = coda.samples(CDF.model, variable.names=c("b00","b01","sigma0","pi0"), 
  #                    numb.iter,thin=numb.iter/5000, progress.bar="none")
  
  
  #########################
  ##### JAGS estimates
  #########################
  b00.hat <- apply(CDF$b00,1:2,mean)
  b01.hat <- apply(CDF$b01,1:2,mean)
  sigma0.hat <- apply(CDF$sigma0,1:2,mean)
  pi0.hat <- apply(CDF$pi0,1:2,mean)
  G <- dim(b00.hat)[2]
  
  b00.hat = apply(b00.hat,1,mean)
  b01.hat = apply(b01.hat,1,mean)
  sigma0.hat = apply(sigma0.hat,1,mean)
  pi0.hat = apply(pi0.hat,1,mean)
  
  ########################
  ###### Estimate PV
  ########################
  
  n1 = length(y1)
  PV.tab <-  rep(NA, n1)
  for(j in 1:n1)
  {
    PV.tab[j] <- 1-sum(pi0.hat*pnorm(y1[j], b00.hat + b01.hat*X1[j], sigma0.hat))
  }
  eta = PV.tab
  
  # inverse = function(fn, interval = NULL, lower = min(interval), upper = max(interval), ...){
  #   Vectorize(function(y){
  #     uniroot(f=function(x){fn(x)-y}, lower=lower, upper=upper, ...)$root
  #   })
  # }
  # F0 = function(c){
  #   sum(pi0.hat*pnorm(c, b00.hat + b01.hat*x, sigma0.hat))
  # }
  # 
  # F0.inv = inverse(F0, lower=-100, upper=100)
  # y1.new = F0.inv(1-eta)
  
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
  
  b0.hat = apply(b0.hat,1,mean)
  b1.hat = apply(b1.hat,1,mean)
  sigma.hat = apply(sigma.hat,1,mean)
  pi.hat = apply(pi.hat,1,mean)
  
  ##################################
  ##### Estimate CDF (ROC) at 
  ##### different iteration level
  ##################################
  grid <- seq(0,1,length=100)

  ROC = sapply(x, function(c){
    
    ROC.tab <-  rep(NA, length(grid))
    for(j in 1:length(grid))
    {
      ROC.tab[j] <- sum(pi.hat*pnorm(qnorm(grid[j]), b0.hat + b1.hat*c, sigma.hat))
    }
    return(ROC.tab)
    })
  AUC = apply(ROC,2,trap)

  ##########################
  ###########
  ##### Optimal functions
  ###########
  ##########################
  
  sp = function(c, x){
    sum(pi0.hat*pnorm(c, b00.hat + b01.hat*x, sigma0.hat))
  }
  
  se = function(c, x){
    sum(pi.hat*pnorm(qnorm(1-sum(pi0.hat*pnorm(c, b00.hat + b01.hat*x, sigma0.hat))), b0.hat + b1.hat*x, sigma.hat))
  }
  
  summary.cutoff = function(c){
    J.c = optim(par = mean(c(y0,y1)), fn = J, method = "BFGS", x = c)
    ER.c = optim(par = mean(c(y0,y1)), fn = ER, method = "BFGS", x = c)
    CZ.c = optim(par = mean(c(y0,y1)), fn = CZ, method = "BFGS", x = c)
    IU.c = optim(par = mean(c(y0,y1)), fn = IU, method = "BFGS", x = c)
    
    J.est = c(J.c$par, -J.c$value, se(J.c$par, x = c), sp(J.c$par, x = c), mean(y1 > J.c$par), mean(y0 < J.c$par))
    ER.est = c(ER.c$par, ER.c$value, se(ER.c$par, x = c), sp(ER.c$par, x = c), mean(y1 > ER.c$par), mean(y0 < ER.c$par))
    CZ.est = c(CZ.c$par, -CZ.c$value, se(CZ.c$par, x = c), sp(CZ.c$par, x = c), mean(y1 > CZ.c$par), mean(y0 < CZ.c$par))
    IU.est = c(IU.c$par, IU.c$value, se(IU.c$par, x = c), sp(IU.c$par, x = c), mean(y1 > IU.c$par), mean(y0 < IU.c$par))
    
    cutoff.tab = rbind(J.est, ER.est, CZ.est, IU.est)
    colnames(cutoff.tab) = c("cutoff","value", "sensitivity", "specificity", 
                             "sensitivity (data)", "specificity (data)")
    return(cutoff.tab)
  }
  
  cutoff = list(
    cutoff0 = summary.cutoff(c = x[1]),
    cutoff1 = summary.cutoff(c = x[2]) #,
    # cutoff3 = summary.cutoff(c = x[3])
  )
    out=list(ROC=ROC, AUC=AUC, cutoff.tab=cutoff)
  return(out)
}


####################
####
####  All fit
####
####################

cutoff.func = function(seed, N = 100, data.gen = BN, x = c(0, 1),
                       burnin = 5000, numb.iter = 5000){
  
  if(N == 50){
    X0 = read.csv("X0_low.csv")$x
    X1 = read.csv("X1_low.csv")$x
  } else if(N == 100){
    X0 = read.csv("X0_medium.csv")$x
    X1 = read.csv("X1_medium.csv")$x
  } else{
    X0 = read.csv("X0_high.csv")$x
    X1 = read.csv("X1_high.csv")$x
  }
  
  dat.mod = do.call(data.gen, list(seed, N, X0, X1, c = x))
  
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
  BN.fit = BiNormal_reg_fit(y0, y1, X0, X1, x, burnin = burnin, numb.iter = numb.iter)
  PV.fit = BayesPV_reg_fit(y0, y1, X0, X1, x, burnin = burnin, numb.iter = numb.iter)
  SemiPV.fit = SemiPar_reg_fit(y0, y1, X0, X1, x, burnin = burnin, numb.iter = numb.iter)
  # Sys.time() - tt
  
  AUC = data.frame(BN = BN.fit$AUC, PV = PV.fit$AUC, Semi.PV = SemiPV.fit$AUC)
  
  BN.ROC = data.frame(grid = grid, BN.fit$ROC) #%>% rename(0 = X1, 1 = X2)
  PV.ROC = data.frame(grid = grid, PV.fit$ROC) #%>% rename(0 = X1, 1 = X2)
  SemiPV.ROC = data.frame(grid = grid, SemiPV.fit$ROC) #%>% rename(0 = X1, 1 = X2)
  ROC = list(BN = BN.ROC, PV = PV.ROC, SemiPV = SemiPV.ROC)

  J.tab = list(cutoff0 = data.frame(BN = BN.fit$cutoff$cutoff0[1,],
                                    PV = PV.fit$cutoff$cutoff0[1,],
                                    Semi.PV = SemiPV.fit$cutoff$cutoff0[1,]),
               cutoff1 = data.frame(BN = BN.fit$cutoff$cutoff1[1,],
                                    PV = PV.fit$cutoff$cutoff1[1,],
                                    Semi.PV = SemiPV.fit$cutoff$cutoff1[1,]))

  ER.tab = list(cutoff0 = data.frame(BN = BN.fit$cutoff$cutoff0[2,],
                                    PV = PV.fit$cutoff$cutoff0[2,],
                                    Semi.PV = SemiPV.fit$cutoff$cutoff0[2,]),
               cutoff1 = data.frame(BN = BN.fit$cutoff$cutoff1[2,],
                                    PV = PV.fit$cutoff$cutoff1[2,],
                                    Semi.PV = SemiPV.fit$cutoff$cutoff1[2,]))
  
  CZ.tab = list(cutoff0 = data.frame(BN = BN.fit$cutoff$cutoff0[3,],
                                     PV = PV.fit$cutoff$cutoff0[3,],
                                     Semi.PV = SemiPV.fit$cutoff$cutoff0[3,]),
                cutoff1 = data.frame(BN = BN.fit$cutoff$cutoff1[3,],
                                     PV = PV.fit$cutoff$cutoff1[3,],
                                     Semi.PV = SemiPV.fit$cutoff$cutoff1[3,]))
  
  IU.tab = list(cutoff0 = data.frame(BN = BN.fit$cutoff$cutoff0[4,],
                                     PV = PV.fit$cutoff$cutoff0[4,],
                                     Semi.PV = SemiPV.fit$cutoff$cutoff0[4,]),
                cutoff1 = data.frame(BN = BN.fit$cutoff$cutoff1[4,],
                                     PV = PV.fit$cutoff$cutoff1[4,],
                                     Semi.PV = SemiPV.fit$cutoff$cutoff1[4,]))
  
  J.tab$cutoff0[1,] = J.tab$cutoff0[1,]*s + mn
  J.tab$cutoff1[1,] = J.tab$cutoff1[1,]*s + mn
  ER.tab$cutoff0[1,] = ER.tab$cutoff0[1,]*s + mn
  ER.tab$cutoff1[1,] = ER.tab$cutoff1[1,]*s + mn
  CZ.tab$cutoff0[1,] = CZ.tab$cutoff0[1,]*s + mn
  CZ.tab$cutoff1[1,] = CZ.tab$cutoff1[1,]*s + mn
  IU.tab$cutoff0[1,] = IU.tab$cutoff0[1,]*s + mn
  IU.tab$cutoff1[1,] = IU.tab$cutoff1[1,]*s + mn
  
  
  out = list(grid = grid, TrueAUC = TrueAUC, TrueROC = TrueROC, True.cutoff = True.cutoff,
             AUC = AUC, ROC = ROC, J.tab = J.tab, ER.tab = ER.tab,
             CZ.tab = CZ.tab, IU.tab = IU.tab)
  return(out)
}


