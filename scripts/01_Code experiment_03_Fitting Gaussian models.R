
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Code Experiments
###' 
###' Task: Fitting the Gaussian Rasch models
###'       
###' Data: Simulated data
###' 
###' Date: 
###' 2022-01-10 created
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

library(brms)
library(edstan)
library(DPpackage)

library(TAM)
library(sirt)
library(lme4)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################
###'
###' Simulate a dataset
###'
###'

nsubject <- 200
nitem <- 40

y <- matrix(0,nrow=nsubject,ncol=nitem)
dimnames(y) <- list(paste("id",seq(1:nsubject)), 
                    paste("item",seq(1,nitem)))


ind <- rbinom(nsubject,1,0.5)
theta <- ind*rnorm(nsubject,1,0.25)+(1-ind)*rnorm(nsubject,3,0.25)
beta <- c(0,seq(-1,3,length=nitem-1))

true.cdf <- function(grid){
  
  0.5*pnorm(grid,1,0.25)+0.5*pnorm(grid,3,0.25) 
  
}  

for(i in 1:nsubject){
  
  for(j in 1:nitem){
    
    eta<-theta[i]-beta[j]         
    
    mean<-exp(eta)/(1+exp(eta))
    
    y[i,j]<-rbinom(1,1,mean)
  }
}


###'#######################################################################
###'
###' Fitting the Bayesian Rasch model with edstan
###'
###'

data_sim <- irt_data(response_matrix = y)

fit_rasch <- irt_stan(data_list = data_sim, iter = 1000, chains = 4)

print_irt_stan(fit_rasch)

as.matrix(fit_rasch)



###'#######################################################################
###'
###' Fitting the Bayesian Rasch model with brms
###'
###'


formula_va_1pl <- bf()








