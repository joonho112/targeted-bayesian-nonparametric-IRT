
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