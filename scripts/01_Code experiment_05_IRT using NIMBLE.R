
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Code Experiments
###' 
###' Task: Item response theory models with NIMBLE
###'       
###' Data: Simulated data
###' 
###' Date: 
###' 2022-03-12
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
work_dir <- file.path(path.expand("~"), 
                      "Documents", 
                      "targeted-bayesian-nonparametric-IRT")

data_dir <- file.path(work_dir, "datasets")

setwd(work_dir)


### Call libraries
library(tidyverse)
library(nimble)
library(igraph)
library(ltm)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################
###'
###' Data
###' 
###' data from the The Law School Admission Test (LSAT) 
###' which is classic example in educational testing (Bock and Lieberman 1970)
###'
###'

LSATdata <- ltm::LSAT

dim(LSATdata)


###'#######################################################################
###'
###' 1PL model
###'
###'

library(nimble, warn.conflicts = FALSE)

code1PL <- nimbleCode({  
  for(i in 1:I) {
    for(p in 1:P) {
      y[p,i] ~ dbern(pi[p,i])
      logit(pi[p, i]) <-  eta[p] - beta[i]
    }
    beta[i] ~ dnorm(0, var = 100)
  }  
  
  for(p in 1:P) {
    eta[p] ~ dnorm(0, sd = sd_eta)
  }
  
  sd_eta ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
})


constants <- list(I = ncol(LSATdata), P = nrow(LSATdata))

data <- list(y = LSATdata)

set.seed(1)
inits <- list(beta    = rnorm(constants$I, 0, 1),
              eta     = rnorm(constants$P, 0, 1), 
              sd_eta = 3)

monitors = c("beta", "eta", "sd_eta")

model1PL <- nimbleModel(code1PL, constants, data, inits)

cModel1PL <- compileNimble(model1PL)


### Building and running an MCMC to fit the model
conf1PL <- configureMCMC(model1PL, monitors = monitors)

model1PLMCMC <- buildMCMC(conf1PL)

cModel1PLMCMC <- compileNimble(model1PLMCMC, project = model1PL)

system.time(samples1PL <- runMCMC(cModel1PLMCMC, niter = 20000, nburnin = 10000))


### Checking the results
betaCols <- grep("beta", colnames(samples1PL))
sd_etaCol <- grep("sd_eta", colnames(samples1PL))
samplesSummary(samples1PL[, c(betaCols, sd_etaCol)])


par(mfrow = c(1, 3), cex = 1.1)
for(i in 1:5)
  ts.plot(samples1PL[ , betaCols[i]], xlab = 'iteration', ylab = colnames(samples1PL)[ betaCols[i]])

ts.plot(samples1PL[ , sd_etaCol], xlab = 'iteration', ylab = colnames(samples1PL)[sd_etaCol])



