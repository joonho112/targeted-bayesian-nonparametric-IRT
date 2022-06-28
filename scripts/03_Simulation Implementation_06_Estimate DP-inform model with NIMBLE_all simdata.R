
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Implementation
###' 
###' Task: Estimate Dirichlet Process Mixture model with informative prior
###'       (DP-inform) using NIMBLE 
###'       across all simulated dataset
###'       
###' Data: Posterior sample data (.rds files)
###' 
###' Date: 
###' - s2022-06-28
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
work_dir <- file.path(path.expand("~"), 
                      "Documents",
                      "targeted-bayesian-nonparametric-IRT") 

data_dir <- file.path(work_dir, "datasets")

data_dir2 <- file.path(path.expand("~"), 
                       "Documents", 
                       "Data-files", 
                       "targeted-bayesian-nonparametric-IRT-large-files")
setwd(work_dir)


### Call libraries
library(tidyverse)

library(nimble)
library(HETOP)

library(future) 
library(furrr)
library(progressr)
library(tictoc)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################
###'
###' Import the pre-generated simulation dataset
###'
###'

### Load the dataset
temp_path <- file.path(data_dir2, 
                       "df_N_items_solutions.rds")

df_pregen <- read_rds(temp_path)


### Subset only one beta generating method
df_pregen %>% 
  count(beta_kind)

df_pre <- df_pregen %>% 
  filter(beta_kind == "beta_norm")

rm(df_pregen)
object.size(df_pre) %>% format("MB")


### Extract an example data
df_example <- df_pre %>%
  filter(
    DGM == "Mixed",
    N_person == 20,
    WLE_rel == 0.5,
    rep == 1
  )

Y <- df_example$y[[1]]
cond_name <- df_example$cond_name[[1]]



###'#######################################################################'
###'
###' Define a function to fit the DPM-inform model & 
###' obtain posterior samples
###'
###'

get_postsamp_DPinform_rasch <- function(Y, cond_name){
  
  # (1) Define a (common) NIMBLE code
  nimble_code <- nimbleCode({ 
    
    for(i in 1:I) {
      for(p in 1:N) {
        y[p,i] ~ dbern(pi[p,i])
        logit(pi[p,i]) <-  theta[p] - beta[i]
      }
    }  
    
    for(i in 1:I) {
      beta[i] ~ dnorm(0,  10)
    } 
    
    # DP mixture model for distribution of ability - CRP representation
    zi[1:N] ~ dCRP(alpha, size = N)
    alpha ~ dgamma(a, b)  
    
    # Mixture component parameter drawn from the base measure
    for(p in 1:N) {
      theta[p] ~ dnorm(mu[p], var = s2[p])  
      mu[p] <- muTilde[zi[p]]                 
      s2[p] <- s2Tilde[zi[p]]   
    }
    
    for(m in 1:M) {
      muTilde[m] ~ dnorm(0, var = s2_mu)
      s2Tilde[m] ~ dinvgamma(nu1, nu2)
    }
    
  })
  
  # (2) Set data, constants, monitors, and inits
  data <- list(y = Y)
  
  constants <- list(I = ncol(Y), N = nrow(Y), M = nrow(Y))
  
  monitors <- c("beta", "theta", 
                "zi", "alpha", 
                # "s2_mu", 
                "a", "b",
                "muTilde", "s2Tilde")
  
  #' Define an informative prior 
  #' with respect to \alpha (the precision parameter)
  #' `N_person = 20, 50, 100, 200, 500`
  info_priors <- list(
    "20" = c(1.24, 0.64),  # df = 5
    "50" = c(1.60, 1.22),  # df = 5
    "75" = c(2.72, 1.36),  # df = 7.5
    "100" = c(3.88, 1.44), # df = 10
    "200" = c(7.80, 1.34), # df = 20
    "300" = c(9.32, 0.88), # df = 30
    "500" = c(9.32, 0.88)  # df = 50
  )
  
  ab_info <- info_priors[[as.character(nrow(Y))]]
  a <- ab_info[1]; b <- ab_info[2]
  
  alpha_mean <- nrow(Y)*0.1
  alpha_var <- a/b^2  
  
  inits <- list(
    beta = rnorm(constants$I, 0, 1), 
    a = a, b = b,
    nu1 = 2.01, nu2 = 1.01, s2_mu = 2, 
    alpha = alpha_mean
  )
  
  scores <- apply(data$y, 1, sum)
  std_scores <- (scores - mean(scores))/sd(scores)
  inits$theta <- std_scores
  inits$zi <- kmeans(std_scores, 5)$cluster
  
  # (3) Create and compile model
  model <- nimbleModel(nimble_code, constants, data, inits)
  
  cModel <- compileNimble(model)
  
  
  # (4) Build an MCMC to fit the model
  conf <- configureMCMC(model, monitors = monitors)
  
  model_MCMC <- buildMCMC(conf)
  
  cModel_MCMC <- compileNimble(model_MCMC, project = model)
  
  
  # (5) Let's roll! - Run an MCMC and obtain posterior samples
  niter <- 4000; nburnin <- 2000
  postsamp <- runMCMC(cModel_MCMC, niter = niter, nburnin = nburnin)
  
  
  #' (6) Save the resulting posterior sample matrix in the designated folder
  #'     data_dir2 need to be present in the global environment
  save_path <-  file.path(data_dir2, 
                          "posterior_sample_DPinform_NIMBLE", 
                          paste0(cond_name, ".rds"))
  
  saveRDS(postsamp, file = save_path)
}


### Define a save version of the function
safe_get_postsamp_DPinform_rasch <- 
  safely(get_postsamp_DPinform_rasch) # safe version

not_null <- negate(is_null)



###'#######################################################################'
###'
###' Fit the DPinform model & save posterior samples 
###' 
###' with parallel computation by future package
###'
###'

### Prepare parallel computation: Set the number of workers
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 10)


### Subset the conditions
df_pre

df_sims_sub <- df_pre

# df_sims_sub <- df_pre %>%
#   filter(N_person %in% c(20, 50))


### Let's roll! - Fit the model
tic()

df_sims_sub %>% 
  # slice(1:10) %>% # slicing for test
  
  ### Fit the models & tidy up results
  future_map2(.x = .$y,
              .y = .$cond_name, 
              .f = ~safe_get_postsamp_DPinform_rasch(.x, .y),
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
  # slice(1:10) %>% # slicing for test
  dplyr::select(sim_cond:beta_kind) %>%
  mutate(error = vec_error)

save_path <- file.path(data_dir2, "df_error_DPinform.rds")

write_rds(df_error, save_path)


### Checking posterior samples invididually
# df_example <- df_sims_sub[10, ]
# 
# (list_est_theta <- check_site_specific_results(postsamp = postsamp,
#                                                param = "theta",
#                                                df_example = df_example))
# 
# (list_est_beta <- check_site_specific_results(postsamp = postsamp,
#                                               param = "beta",
#                                               df_example = df_example))
# 
# 
# ### DPM model check
# DPdiffuse_check <- DPM_model_check(postsamp = postsamp)
# 
# DPdiffuse_check[1]
# DPdiffuse_check[2]
# DPdiffuse_check[3]

