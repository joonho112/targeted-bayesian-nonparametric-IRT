###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Code Experiments
###' 
###' Task: Simulating data from DPpackage & sirt R package 
###'       
###' Data: Simulated data
###' 
###' Date: 2021-10-16
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


### Set working directory and data directory 
work_dir <- c("~/Documents/targeted-bayesian-nonparametric-IRT")
data_dir <- file.path(work_dir, "datasets")
setwd(work_dir)


### Call libraries
library(tidyverse)
library(TAM)
library(sirt)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'######################################################################
###'
###' Simulated data example from `DPpackage`
###'
###'

nsubject <- 250

nitem <- 40

y <- matrix(0, nrow = nsubject, ncol = nitem) 

dimnames(y)<-list(
  paste("id",seq(1:nsubject)), 
  paste("item",seq(1, nitem)))

ind <- rbinom(nsubject, 1, 0.5) 

theta <- ind*rnorm(nsubject, -1, sqrt(0.25)) + 
  (1 - ind)*rnorm(nsubject, 2, sqrt(0.065)) 

beta <- c(0, seq(-3, 3, length = nitem-1))

true.density <- function(grid) { 
  0.5*dnorm(grid, -1, sqrt(0.25)) + 0.5*dnorm(grid, 2, sqrt(0.065)) 
  }

true.cdf <- function(grid) { 
  0.5*pnorm(grid, -1, sqrt(0.25)) +0.5*pnorm(grid, 2, sqrt(0.065)) 
  }

for(i in 1:nsubject) { 
  for(j in 1:nitem) { 
    eta <- theta[i] - beta[j] 
    prob <- exp(eta)/(1 + exp(eta)) 
    y[i,j] <- rbinom(1, 1, prob) 
  } 
}



###'######################################################################
###'
###' Simulated data example from `sirt::sim.raschtype()`
###'
###' - Simulation of data from a Rasch model (alpha_1=alpha_2=0)
###'
###'

### Simulate setup
set.seed(9875)
N <- 500  # number of persons
I <- 11   # number of items
b <- sample(seq(-2, 2, length = I))
a <- rep(1, I)


### Create some misfitting items 
a[c(1, 3)] <- c(.5, 1.5)


### Simulate data
dat <- sirt::sim.raschtype(theta = rnorm(N), b = b, fixed.a = a)



###'######################################################################
###'
###' Extracting individual likelihood and individual posterior
###'
###' - Dichotomous data data.sim.rasch
###'
###'

### Load the data
data("data.sim.rasch")


### 1PL estimation
mod1 <- TAM::tam.mml(resp = data.sim.rasch)


### Extract likelihood
lmod1 <- IRT.likelihood(mod1)

str(lmod1)


### Extract posterior
pmod1 <- IRT.posterior(mod1)

str(pmod1)

df_temp <- pmod1 %>% data.frame()

dim(df_temp)



###'#######################################################################'
###'
###' Toy example for reliability functions
###'
###'

set.seed(9897)

N <- 100

### Simulate theta and error SDs
x <- stats::rnorm(N, sd = 2)

error <- stats::runif(N, 0.7, 1.3)


### Compute WLE reliability
WLErel(x, error)


### Compute EAP reliability
EAPrel(x, error)


