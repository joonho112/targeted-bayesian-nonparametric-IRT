
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Data generation 
###' 
###' Task: Examining the simulated data
###'       generated according to the level of reliabilities
###'
###' Data: Simulated data
###' 
###' Date: 
###' - 2021-12-06
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
data_dir2 <- c("~/Documents/Data-files/targeted-bayesian-nonparametric-IRT-large-files")
setwd(work_dir)


### Call libraries
library(tidyverse)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################'
###'
###' Load the simulated data 
###'
###'

setwd(data_dir2)

df <- readRDS(file = "simulation_results.rds")


### Check up the data structure
df$DGM

df$N_person

df$WLE_rel

df$true_theta[[1]]

df$grid_search[[1]]

df$solution[[1]]

df$solution_sub[[1]]


### Factor levels 
vec_DGM <- c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")


###'#######################################################################'
###'
###' The number of items (`I`) per each `WLE` reliability
###' 
###'

setwd(work_dir)

###' Unnest the solution_sub tibble
df_sub <- df %>%
  select(solution_sub) %>%
  unnest(solution_sub) %>%
  mutate(DGM = factor(DGM, levels = vec_DGM))


###' The average number of items
###' regardless of item difficulty generation functions
###' 
###' => `N_person` and `DGM` doesn't influence the average I level
###' => Only `WLE reliability` and `beta_kind` matter
###' => This is because we marginalized over the person dimension
###'    to obtain item information. 
df_I_solution <- df_sub %>%
  group_by(DGM, N_person, WLE_rel) %>%
  summarize(I_mean = mean(I)) %>%
  arrange(WLE_rel, N_person, DGM) %>%
  relocate(DGM, .after = WLE_rel)


###' The average number of items
###' accroding to item difficulty generation functions
###' => Again, I varies only across `WLE_reliability` and `beta_kind`
df_I_solution_beta <- df_sub %>%
  group_by(DGM, beta_kind, N_person, WLE_rel) %>%
  summarize(I_mean = mean(I)) %>%
  arrange(WLE_rel, N_person, DGM, beta_kind) %>%
  relocate(DGM, beta_kind, .after = WLE_rel)


###' Thus, now we can collapse the dimensions of `N_person` and `DGM`

df_I_solution_beta2 <- df_sub %>%
  group_by(WLE_rel, beta_kind) %>%
  summarize(I_mean = mean(I)) %>%
  arrange(WLE_rel, beta_kind)

df_I_solution_beta3 <- df_I_solution_beta2 %>%
  pivot_wider(names_from = beta_kind, values_from = I_mean) # Final solution


