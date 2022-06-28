
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Implementation
###' 
###' Task: Build a funcion to 
###'       fit the Dirichlet Process Mixture Rasch models
###'       `with diffuse prior`
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2022-02-18 finalized
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


### Set working directory and data directory (Mac) 
work_dir <- c("~/Documents/targeted-bayesian-nonparametric-IRT")
data_dir <- file.path(work_dir, "datasets")
data_dir2 <- c("~/Documents/Data-files/targeted-bayesian-nonparametric-IRT-large-files")
setwd(work_dir)


### Set working directory and data directory (Windows)
work_dir <- c("~/targeted-bayesian-nonparametric-IRT")
data_dir <- file.path(work_dir, "datasets")
data_dir2 <- c("D:/Data-files/targeted-bayesian-nonparametric-IRT-large-files")
setwd(work_dir)


### Call libraries
library(tidyverse)

library(DPpackage)
library(TAM)
library(sirt)

library(future) 
library(furrr)
library(progressr)
library(tictoc)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'######################################################################
###'
###' Load the pre-generated simulation data
###' 
###' - `df_N_items_solutions`: 
###'
###'

### Load the dataset
temp_path <- file.path(data_dir2, 
                       "df_N_items_solutions.rds")

df <- read_rds(temp_path)


### Subset only one beta generating method
df %>% 
  count(beta_kind)

df_sub <- df %>% 
  filter(beta_kind == "beta_norm")


# ### Extract the example data
# y_mat <- df$y[[1000]]



###'#######################################################################'
###'
###' Define a function to fit the DP-diffuse model & 
###' obtain posterior samples
###'
###'

get_postsamp_DPdiffuse_rasch <- function(y_mat, filename){
  
  # STEP 1: Set MCMC parameters
  N_iter <- 4000
  state <- NULL    # the current value of the parameters
  mcmc <- list(nburn = N_iter,      # the number of burn-in scans, 
               nsave = N_iter,      # the total number of scans to be saved,
               nskip = N_iter/200,  # the thinning interval, 
               ndisplay = N_iter/40)   # number of saved scans displayed on screen
  
  
  # STEP 2: Define a Diffuse prior 
  # with respect to \alpha (the precision parameter)
  n_person <- nrow(y_mat)
  n_item <- ncol(y_mat)

  # (1) alpha precision parameter
  b <- 0.1
  alpha_mean <- n_person/2
  a <- alpha_mean*b
  alpha_var <- a/b^2  
  
  # (2) beta parameters
  beta0 <- rep(0, n_item - 1)
  Sbeta0 <- diag(100, n_item - 1)

  prior <- list(a0 = a,      # alpha0: shape param.
                b0 = b,      # alpha0: rate param.  
                tau1 = 1,    # G0 variance: shape param.
                tau2 = 1,    # G0 variance: rate param. 
                mub = 0,     # G0 mean: mean param. 
                Sb = 100,    # G0 mean: variance param. 
                beta0 = beta0, 
                Sbeta0 = Sbeta0)    # G0 mean: mean 

  
  # STEP 3: Fit the Rasch model with the DP diffuse prior
  outp <- DPrasch(y = y_mat, 
                  prior = prior, 
                  mcmc = mcmc, 
                  state = state, 
                  status = TRUE)
  
  
  # STEP 4: Extract posterior samples
  list_posterior <- get_posterior_DPrasch(output = outp, nburn_DP = N_iter)
  
  
  # STEP 5: Save the resulting list in the designated folder
  save_path <-  file.path(data_dir, 
                          paste0(filename, ".rds"))
  
  saveRDS(list_posterior, file = save_path)
  
  # return(df_posterior)
}

safe_postsamp_DPdiffuse_rasch <- 
  safely(get_postsamp_DPdiffuse_rasch) # safe version

not_null <- negate(is_null)




###'#######################################################################'
###'
###' Fit the Gaussian model & save posterior samples 
###' 
###' with parallel computation by future package
###'
###'

### Prepare parallel computation: Set the number of workers
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 20)


### Subset the conditions
df_sub

df_sims_sub <- df_sub %>%
  filter(N_person %in% c(100))


### Let's roll! - Fit the models
tic()

df_sims_sub %>% 
  #slice(1:10) %>% # slicing for test
  
  ### Fit the models & tidy up results
  future_map2(.x = .$y,
              .y = .$cond_name, 
              .f = ~safe_postsamp_DPdiffuse_rasch(.x, .y),
              .options = furrr_options(seed = NULL,
                                       chunk_size = NULL,
                                       scheduling = 10),
              .progress = TRUE) %>%
  map("error") %>%
  map_lgl(not_null) %>%
  {. ->> vec_error}

toc()


### Check out errors and save the dataframe
df_error <- df_sims_sub %>%
  #slice(1:10) %>% # slicing for test
  dplyr::select(sim_cond:beta_kind) %>%
  mutate(error = vec_error)

save_path <- file.path(data_dir, "df_error_DP-diffuse_100.rds")

write_rds(df_error, save_path)
  






