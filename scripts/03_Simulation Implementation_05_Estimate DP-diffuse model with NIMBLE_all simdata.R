
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Implementation
###' 
###' Task: Estimate Dirichlet Process Mixture model with diffuse prior
###'       (DP-diffuse) using NIMBLE 
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
                      # "Documents",
                      "targeted-bayesian-nonparametric-IRT") 

data_dir <- file.path(work_dir, "datasets")

data_dir2 <- file.path(path.expand("~"), 
                       # "Documents", 
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
###' Define a function to fit the DPM-diffuse model & 
###' obtain posterior samples
###'
###'

get_postsamp_DPdiffuse_rasch <- function(Y, cond_name){
  
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
  
  b <- 0.1
  alpha_mean <- nrow(Y)/2
  a <- alpha_mean*b
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
                          "posterior_sample_DPdiffuse_NIMBLE", 
                          paste0(cond_name, ".rds"))
  
  saveRDS(postsamp, file = save_path)
}


### Define a save version of the function
safe_get_postsamp_DPdiffuse_rasch <- 
  safely(get_postsamp_DPdiffuse_rasch) # safe version

not_null <- negate(is_null)



###'#######################################################################'
###'
###' Fit the DPdiffuse model & save posterior samples 
###' 
###' with parallel computation by future package
###'
###'

### Prepare chopped looping for memory savings
vec_DGM <- c("Gaussian", "ALD", "Mixed")
vec_N_person <- c(20, 50, 100, 200, 500)
<<<<<<< HEAD
# list_rep <- split(51:100, rep(1:17, each = 3))
list_rep <- split(70:81, rep(1:12, each = 1))
=======
list_rep <- split(1:81, rep(1:27, each = 3))
>>>>>>> e90deac1a6bf28fb11481ad80f9bf2b31a119efd


### Construct a for loop for memory savings
tic()

for (i in seq_along(list_rep)){
  for (j in seq_along(vec_DGM)){
    for (k in seq_along(vec_N_person)){
      
      # Print progress
      cat(paste0("DGM: ", vec_DGM[j]), 
          paste0("N_person:", vec_N_person[k]), "rep: ", 
          paste0(list_rep[[i]]), "\n")
      
      # Prepare parallel computation: Set the number of workers
      # parallelly::availableCores()
      # parallelly::availableWorkers()
<<<<<<< HEAD
      plan(multisession, workers = 5)
=======
      plan(multisession, workers = 15)
>>>>>>> e90deac1a6bf28fb11481ad80f9bf2b31a119efd
      
      # Subset the conditions
      df_sims_sub <- df_pre %>%
        filter(DGM %in% vec_DGM[j]) %>%
        filter(N_person %in% vec_N_person[k]) %>%
        filter(rep %in% list_rep[[i]])
      
      # Let's roll! - Fit the model
      tic()
      
      df_sims_sub %>% 
        # slice(1:2) %>% # slicing for test
        
        ### Fit the models & tidy up results
        future_map2(.x = .$y,
                    .y = .$cond_name, 
                    .f = ~safe_get_postsamp_DPdiffuse_rasch(.x, .y),
                    .options = furrr_options(seed = NULL,
                                             chunk_size = NULL,
                                             scheduling = 2),
                    .progress = TRUE) %>%
        map("error") %>%
        map_lgl(not_null) %>%
        {. ->> vec_error}
      
      toc()
      
    }
  }
}

toc()


<<<<<<< HEAD

=======
>>>>>>> e90deac1a6bf28fb11481ad80f9bf2b31a119efd

# ### Checking posterior samples invididually
# df_example <- df_sims_sub[8, ]
# 
# (list_est_theta <- check_site_specific_results(postsamp = postsamp,
#                                                param = "theta",
#                                                df_example = df_example))
# 
# (list_est_beta <- check_site_specific_results(postsamp = postsamp,
#                                               param = "beta",
#                                               df_example = df_example))


# ### DPM model check
# DPdiffuse_check <- DPM_model_check(postsamp = postsamp)
# 
# DPdiffuse_check[1]
# DPdiffuse_check[2]
# DPdiffuse_check[3]

