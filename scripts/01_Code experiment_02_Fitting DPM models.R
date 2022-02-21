
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Code Experiments
###' 
###' Task: Fitting the Dirichlet Process Mixture Rasch models
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2022-01-10 created
###' - 2022-02-18 revised
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###' 

###'######################################################################
###'
###' Basic settings
###'
###'

### Start with a clean slate (release memory; delete objects)
gc(); rm(list=ls())


### Set working directory and data directory (for Mac)
work_dir <- c("~/Documents/targeted-bayesian-nonparametric-IRT")
data_dir <- file.path(work_dir, "datasets")
setwd(work_dir)


### Set working directory and data directory (for Windows)
work_dir <- c("~/targeted-bayesian-nonparametric-IRT")
data_dir <- file.path(work_dir, "datasets")
setwd(work_dir)


### Call libraries
library(tidyverse)

library(DPpackage)
library(TAM)
library(sirt)

library(tictoc)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################
###'
###' `DPrasch()`
###' - Bayesian analysis for a semiparametric Rasch model
###' 
###'

### A simulated Data Set
nsubject <- 200
nitem <- 40

y <- matrix(0, nrow = nsubject, ncol = nitem)

dimnames(y) <- list(paste("id", seq(1:nsubject)), 
                    paste("item", seq(1,nitem)))


ind <- rbinom(nsubject, 1, 0.5)

theta <- ind*rnorm(nsubject, 1, 0.25) + (1-ind)*rnorm(nsubject, 3, 0.25)

beta <- c(0, seq(-1, 3, length = nitem - 1))

true.cdf <- function(grid){
  
  0.5*pnorm(grid,1,0.25) + 0.5*pnorm(grid,3,0.25) 

}  

for(i in 1:nsubject){
  
  for(j in 1:nitem){
    
    eta <- theta[i]-beta[j]         
    
    mean <- exp(eta)/(1+exp(eta))
    
    y[i,j] <- rbinom(1, 1, mean)
  }
}


### Prior information
beta0 <- rep(0, nitem - 1)
Sbeta0 <- diag(1000, nitem - 1)

prior <- list(alpha = 1,  # instead of a0, b0
              tau1 = 6.02,
              tau2 = 2.02,
              mub = 0,
              Sb = 100,
              beta0 = beta0,
              Sbeta0 = Sbeta0)

### Initial state
state <- NULL


### MCMC parameters
nburn <- 2000
nsave <- 2000
nskip <- 0
ndisplay <- 1000

mcmc <- list(nburn = nburn,
             nsave = nsave,
             nskip = nskip,
             ndisplay = ndisplay)


### Fit the model
tic()

fit1 <- DPrasch(y = y,
                prior = prior, 
                mcmc = mcmc, 
                state = state,
                status = TRUE,
                grid = seq(-1, 5, 0.01),
                compute.band = TRUE)

toc()


### CDF estimate and truth
plot(fit1$grid, 
     true.cdf(fit1$grid),
     type="l", lwd = 2, col = "red",
     xlab = expression(theta),
     ylab = "CDF")

lines(fit1$grid, fit1$cdf, 
      lwd = 2, col = "blue")

lines(fit1$grid, fit1$cdf.l, 
      lwd = 2, col = "blue", lty = 2)

lines(fit1$grid, fit1$cdf.u, 
      lwd = 2, col = "blue", lty = 2)


### Summary with HPD and Credibility intervals
summary(fit1)
summary(fit1, hpd = FALSE)


###' Plot model parameters 
###' (to see the plots gradually set ask=TRUE)
plot(fit1, ask = FALSE)
plot(fit1, ask = FALSE, nfigr = 2, nfigc = 2)	


###' Extract random effects
DPrandom(fit1)
plot(DPrandom(fit1))
DPcaterpillar(DPrandom(fit1))


### Extract posterior samples
# (1) Tidy up the "theta" output object
df_theta <- fit1$save.state$randsave %>%
  data.frame() %>% 
  dplyr::select(-theta.Prediction.) %>%
  slice(c((nburn/2 + 1):nburn)) %>%  # throw away burn-ins
  set_names(paste0("theta[", seq(ncol(.)), "]")) %>%
  tibble()


# (2) Tidy up the G0 parameters (m, s2), alpha0, and N of clusters outputs
df_hyperparm <- fit1$save.state$thetasave %>%
  data.frame() %>%
  dplyr::select(-tau_j_hat) %>%
  rename(G0_mu = mu, G0_s2 = sigma2, Ncluster = ncluster, alpha0 = alpha) %>%
  dplyr::select(G0_mu, G0_s2, alpha0, Ncluster) %>%
  slice(c((nburn_DP/2 + 1):nburn_DP))  # throw away burn-ins

# Return the resulting object
cbind.data.frame(df_hyperparm, df_theta)



###'#######################################################################
###'
###' `DPMrasch()`
###' 
###'  Bayesian analysis for a semiparametric Rasch model
###'
###'

### A simulated Data Set
nsubject <- 250
nitem <- 40

y <- matrix(0, nrow = nsubject, ncol = nitem)

dimnames(y)<-list(paste("id",seq(1:nsubject)), 
                  paste("item",seq(1,nitem)))

ind <- rbinom(nsubject, 1, 0.5)

theta <- ind*rnorm(nsubject, -1, sqrt(0.25))+
  (1-ind)*rnorm(nsubject ,2, sqrt(0.065))

beta <- c(0, seq(-3, 3, length = nitem - 1))


true.density <- function(grid){
  0.5*dnorm(grid, -1, sqrt(0.25)) + 0.5*dnorm(grid, 2, sqrt(0.065))  
} 

true.cdf <- function(grid){
  0.5*pnorm(grid, -1, sqrt(0.25)) + 0.5*pnorm(grid, 2, sqrt(0.065))  
} 

for(i in 1:nsubject)
{
  for(j in 1:nitem)
  {
    eta <- theta[i] - beta[j]         
    prob <- exp(eta)/(1 + exp(eta))
    y[i,j] <- rbinom(1, 1, prob)
  }
}

### Prior information
beta0 <- rep(0, nitem - 1)
Sbeta0 <- diag(100, nitem - 1)

prior <- list(N = 50, # ? not clear
              alpha = 1,
              taub1 = 6.01,
              taub2 = 2.01,
              taus1 = 6.01,
              taus2 = 2.01,
              tauk1 = 6.01,
              m0= 0,
              s0 = 100,
              beta0 = beta0,
              Sbeta0 = Sbeta0)

### Initial state
state <- NULL      


### MCMC parameters
nburn <- 4000
nsave <- 4000
nskip <- 0
ndisplay <- 100
mcmc <- list(nburn = nburn,
             nsave = nsave,
             nskip = nskip,
             ndisplay = ndisplay)


### Fit the model
tic()

fit1 <- DPMrasch(y = y,
                 prior = prior,
                 mcmc = mcmc,
                 state = state,
                 status = TRUE,
                 grid = seq(-3, 4, 0.01))

toc()


### Plot the results
plot(fit1$grid, fit1$dens.m,
     type="l", lty = 1, col = "red",
     xlim = c(-3,4), 
     ylim = c(0,0.8))

lines(fit1$grid,true.density(fit1$grid),
      lty = 2, col = "blue")

plot(fit1$grid, 
     fit1$cdf.m, 
     type = "l", lty = 1, col = "red")

lines(fit1$grid, 
      true.cdf(fit1$grid),
      lty = 2, col = "blue")


### Summary with HPD and Credibility intervals
summary(fit1)
summary(fit1,hpd=FALSE)


###' Plot model parameters 
###' (to see the plots gradually set ask=TRUE)
plot(fit1, ask = FALSE)
plot(fit1, ask = FALSE, nfigr = 2, nfigc = 2)	


### Extract random effects
DPrandom(fit1)
plot(DPrandom(fit1))
DPcaterpillar(DPrandom(fit1))

