###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Code Experiments
###' 
###' Task: Simulating data from Dr. Stefanie Wind's code
###'       
###' Data: Simulated data
###' 
###' Date: 2021-10-25
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



GPR<- function(N=NULL,model="GPCM",theta=NULL,discrim=NULL,delta=NULL,tau=NULL){
  m2l <- function(m,remove=NA){
    if(is.na(remove)){
      lapply(seq_len(nrow(m)), function(i) m[i,!is.na(m[i, ]) ])
    }else{
      lapply(seq_len(nrow(m)), function(i) m[i,m[i, ]!=remove ])
    }
  }
  J=length(delta)
  if(is.null(theta)){
    theta=rnorm(N)
  }else if(is.null(N)){
    N=length(theta)
  }
  if(model=="PCM"){
    discrim <- rep(1,J)
  }else if(model=="RSM"){
    discrim <- rep(1,J)
    delta <- outer(delta,tau,"+")
    delta<- m2l(delta)
  } 
  resp_all <- NULL
  for(i in c(1:J)){ ##i is number of items
    theta_d <- discrim[i]*cbind(0,outer(theta,delta[[i]],"-"))
    cumsum_t_d <- t(apply(theta_d,1,cumsum))
    exp1 <- exp(cumsum_t_d)
    p <- exp1/rowSums(exp1)
    resp0 <- NULL
    for(n in c(1:N)){## n is the number of students
      resp <- sample(0:length(delta[[i]]),1,prob =p[n,])## 1st student's response to item 1 
      resp0<- c(resp0,resp)
    }
    resp_all<- cbind(resp_all,resp0)
    colnames(resp_all) <- NULL
  }
  return(resp_all)
}


###Examples###
###use GPR function to generate data based on RSM###
N=2000 ## sample size = 2000, if users don't provide sample size, they need to provide theta parameter. It is OK if users would like to provide both sample size and theta parameter.
model = "RSM"
delta=c(-1,-1.2,0,1,1.5,2) ## 6 items
tau=c(-1,0,1) ## each item has 4 categories. Only in RSM, users should provide tau parameter.
data <- GPR(N,model="RSM",delta=delta,tau=tau)

###use GPR function to generate data based on PCM###
model = "PCM"
theta=rnorm(1500,0,1.5)
delta=list(c(-1,0),c(-1,1),c(-1,0,1),c(-2,-0.5,0.5,1.3),c(1,0),c(1.5,-0.5)) ## 6 items, each item can have different the number of categories. The second step is not necessary harder than the first step.
data <- GPR(model="PCM",theta=theta,delta=delta)

###use GPR function to generate data based on GPCM###
N=2000
discrim=runif(6,0.4,2) ## 6 items, each item has different discrimination parameters. Users can also provide the same discrimination parameter for each item.
delta=list(c(-1,0),c(-1,1),c(-1,0,1),c(-2,-0.5,0.5,1.3),c(1,0),c(1.5,-0.5))
data <- GPR(N,discrim=discrim,delta=delta)

