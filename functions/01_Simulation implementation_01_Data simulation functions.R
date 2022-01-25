
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Define functions
###' 
###' Task: Define required functions
###'       to simulate item response data according to
###'       the level of reliabilities
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2021-11-14 
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###' 

###'#######################################################################'
###'
###' Call necessary packages
###'
###'

library(tidyverse)
library(kableExtra)
library(psych)
library(sirt)
library(TAM)
library(LaplacesDemon)
library(sn)



###'#######################################################################'
###' 
###' `gen_sim_conds()`
###' 
###'  - generate a table with simulation conditions
###'    depending upon true variance and WLE reliability
###'
###'

gen_sim_conds <- function(true_var, WLE_rel){
  
  # True SD
  true_SD <- sqrt(true_var)
  
  # Observed variance = true variance + error variance
  obs_var <- true_var/WLE_rel
  
  # Error variance and RMSE
  err_var <- obs_var - true_var
  err_SD <- sqrt(err_var)
  
  # Signal-to-noise ratio = Separation coefficient = true SD / error SD (RMSE)
  sep_coef <- true_SD/err_SD
  
  # Strata = (4*sep_coef + 1)/3
  strata <- (4*sep_coef + 1)/3
  
  # Test information = 1/ error variance
  test_info <- 1/err_var
  
  # EAP reliability
  EAP_rel <- 1 - (err_var/(obs_var + err_var))
  
  # Combine as a tibble
  tab_cond <- tibble(
    true_SD, true_var, err_SD, err_var, obs_var, 
    sep_coef, strata, test_info, WLE_rel, EAP_rel
  )
  
  return(tab_cond)
}


# ### Testing
# 
# # Target WLE reliabilities and fixed true variances
# WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# true_var <- c(1, 1, 1, 1, 1)
# 
# # New simulation conditions
# df_simconds <- gen_sim_conds(true_var, WLE_rel)
# 
# round(df_simconds, 2) %>%
#   kbl(caption = "Simulation conditions #02") %>%
#   kable_minimal(full_width = F, html_font = "Cambria")



###'#######################################################################'
###'
###' `tibble_sim_conditions()`
###' 
###' Simulation tibble generator
###' 
###' - Number of persons (test takers)
###' - Data-generating functions for the latent trait (\theta) distribution
###' - true variance of \theta
###' - WLE reliability
###' 
###'

tibble_sim_conditions <- function(vec_N_persons, 
                                  vec_DGMs, 
                                  vec_true_var, 
                                  vec_WLE_rel){
  
  # Generate simulation conditions for WLE reliability
  df_simconds <- gen_sim_conds(vec_true_var, vec_WLE_rel)
  
  # A table for all combinations - person-side dimension 
  df_theta_conds <- list(
    DGM = vec_DGMs, 
    N_person = vec_N_persons
  ) %>% cross_df()
  
  # Combine WLE reliability table and person-side conditions
  df_theta_conds <- df_theta_conds %>%
    mutate(sim_cond = rep(list(df_simconds), nrow(df_theta_conds))) %>%
    unnest(sim_cond)
  
  return(df_theta_conds)
} 


# ### Test the function
# vec_N_persons <- c(20, 50, 100, 200, 1000)
# vec_DGMs <- c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")
# vec_WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# vec_true_var <- c(1, 1, 1, 1, 1)
# 
# df_simcond <- tibble_sim_conditions(vec_N_persons, vec_DGMs, vec_true_var, vec_WLE_rel)



###'#######################################################################'
###'
###' `gen_true_theta()`
###' 
###'  - generate theta's: true latent trait distributions
###'
###'

### Define a function to generate theta distribution
gen_true_theta <- function(N = 50, true_dist = NULL, 
                           tau = 0, var = 1, 
                           nu = NULL, slant = NULL, rho = NULL, 
                           delta = NULL, eps = NULL, ups = NULL){
  
  # Set mean, variance, and SD 
  # zero mean and unit variance for all Gs   
  sigma <- sqrt(var)
  
  # Generate true person ability distribution theta
  
  if (true_dist == "Gaussian"){ 
    
    theta <- rnorm(N, mean = tau, sd = sigma)
    
  } else if (true_dist == "T" & !is.null(nu)){   
    
    theta <- rt(N, nu)*sqrt((nu - 2)/nu)
    
  } else if (true_dist == "Skew" & !is.null(slant)){
    
    # Generate location and scale parameters of skewed normal
    # to achieve E(theta) = 0 and Var(theta) = 1
    delta <- slant/sqrt(1 + slant^2)
    scale <- sqrt(var/(1-(2*delta^2/pi)))
    location <- tau - scale*sqrt(2/pi)*delta
    
    theta <- sn::rsn(n = N, xi = location, omega = scale, 
                     alpha = slant)[1:N]
    
  } else if (true_dist == "ALD" & !is.null(rho)){
    
    # Generate location, scale, and skewness parameters 
    # of ALD distribution to achieve E(theta) = 0 and Var(theta) = 1
    scale <- sqrt((2*rho^2*var)/(1 + rho^4))
    location <- tau - ((scale*(1/rho - rho))/sqrt(2))
    theta <- LaplacesDemon::ralaplace(N, location, scale, rho)
    
  } else if (true_dist == "Mixture" & !(is.null(delta)|is.null(eps)|is.null(ups))){
    
    # Define a normalizing factor `a`
    a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
    
    # Simulate a mixture of two normals with mean 0 and variance 1  
    ind <- runif(N) < (1 - eps)
    theta <- ind*rnorm(N, -eps*delta/a, sqrt(1/a^2)) + 
      (1 - ind)*rnorm(N, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  }
  return(theta)
}

# ### Testing
#
# ### Simulate theta's according to the defined DGMs
# N <- 10000
# 
# Gaussian <- gen_true_theta(N = N, true_dist = "Gaussian", tau = 0, var = 1)
# 
# T <- gen_true_theta(N = N, true_dist = "T", tau = 0, var = 1, nu = 5)
# 
# Skew <- gen_true_theta(N = N, true_dist = "Skew", tau = 0, var = 1, slant = 5)
# 
# ALD <-gen_true_theta(N = N, true_dist = "ALD", tau = 0, var = 1, rho = 0.1)
# 
# Bimodal <- gen_true_theta(N = N, true_dist = "Mixture", tau = 0, var = 1,
#                           delta = 4, eps = 0.3, ups = 1)
# 
# Mixed <- gen_true_theta(N = N, true_dist = "Mixture", tau = 0, var = 1,
#                         delta = 5, eps = 0.3, ups = 2)
# 
# ### Collect the simulated data
# vec_levels <- c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")
# 
# df_theta <- tibble(Gaussian, T, Skew, ALD, Bimodal, Mixed) %>%
#   pivot_longer(Gaussian:Mixed, names_to = "DGM", values_to = "value") %>%
#   mutate(DGM = factor(DGM, levels = vec_levels)) %>%
#   arrange(DGM)
# 
# df_theta



###'#######################################################################'
###'
###' `gen_vec_theta()`
###'
###' - a wrapper function to generate theta according to DGM and N_person
###'
###'

gen_vec_theta <- function(DGM, N){
  
  if (DGM == "Gaussian"){ 
    vec_theta <- gen_true_theta(N = N, true_dist = "Gaussian", tau = 0, var = 1)
    
  } else if (DGM == "T"){   
    vec_theta <- gen_true_theta(N = N, true_dist = "T", tau = 0, var = 1, 
                                nu = 5)
    
  } else if (DGM == "Skew"){
    vec_theta <- gen_true_theta(N = N, true_dist = "Skew", tau = 0, var = 1, 
                                slant = 5)
    
  } else if (DGM == "ALD"){
    vec_theta <- gen_true_theta(N = N, true_dist = "ALD", tau = 0, var = 1, 
                                rho = 0.1)
    
  } else if (DGM == "Bimodal"){
    vec_theta <- gen_true_theta(N = N, true_dist = "Mixture", tau = 0, var = 1,
                                delta = 4, eps = 0.3, ups = 1)
    
  } else if (DGM == "Mixed"){
    vec_theta <- gen_true_theta(N = N, true_dist = "Mixture", tau = 0, var = 1,
                                delta = 5, eps = 0.3, ups = 2)
  }
}



###'#######################################################################'
###'
###' `gen_rasch_logit_prob_y()`
###'
###' - generate logit, probability, and binary responses 
###' 
###'

gen_rasch_logit_prob_y <- function(N_person, N_item, theta, beta){
  
  # Create two empty matrices saving probabilities and responses 
  y_mat <- prob_mat <- logit_mat <- matrix(0, nrow = N_person, ncol = N_item) 
  
  # Generate responses based on Rasch model
  for(i in 1:N_person) { 
    for(j in 1:N_item) { 
      
      logit <- theta[i] - beta[j]
      logit_mat[i, j] <- logit
      
      prob <- exp(logit)/(1 + exp(logit))
      prob_mat[i, j] <- prob
      
      y_mat[i,j] <- rbinom(1, 1, prob) 
    } 
  }
  
  # Return three matrices as a list
  dimnames(y_mat) <- dimnames(prob_mat) <- dimnames(logit_mat) <-
    list(
      paste0("id_", str_pad(seq(1:N_person), 2, pad = "0")), 
      paste0("item_", str_pad(seq(1:N_item), 2, pad = "0"))
    )
  list(logit_mat, prob_mat, y_mat)
}

# ### Testing
# N <- 100
# I <- 20
# 
# theta <- gen_true_theta(N = N, true_dist = "Gaussian", tau = 0, var = 1)
# beta <- c(0, seq(-2, 2, length = I - 1)) %>% sort()
# 
# list_logit_prob_y <- gen_rasch_prob_y(N, I, theta, beta)
# 
# list_logit_prob_y



###'#######################################################################
###'
###' `get_test_info_est()`
###' 
###'  - Calculate test information from the logit matrix
###'
###'

get_test_info_est <- function(logit_mat){
  
  # Item information vector
  logit_item <- apply(logit_mat, 2, mean)
  
  # Test information scalar (estimated)
  prob_item <- exp(logit_item)/(1 + exp(logit_item))
  item_info <- prob_item * (1 - prob_item)
  test_info_est <- sum(item_info)
  return(test_info_est)
}



###'#######################################################################'
###'
###' `prophecy_target_N_item()`
###'
###'  - Generate the target number of items
###'    based on the Spearman-Brown prediction formula 
###'    for person test reliability
###'
###'

prophecy_target_N_item <- function(C = 20, R_C = 0.88, R_T = 0.50){
  
  T <- C*(R_T*(1 - R_C))/((1 - R_T)*R_C)
  return(T)
  
}


###'#######################################################################'
###'
###' `fit_rasch()`
###'
###' - Fit Rasch model 
###'
###'

fit_rasch <- function(y_mat){
  
  # Fit a simple Rasch model (MML estimation) 
  TAM::tam.mml(resp = y_mat)
  
}
  

###'#######################################################################'
###'
###' `get_reliabilities()`
###'
###' - Calculate WLE and EAP reliabilities from the simulated binary responses  
###' 
###'
 
get_reliabilities <- function(mod){
  
  # WLE estimation
  wle <- TAM::tam.wle(mod)
  
  # Return WLE and EAP reliabilities
  tibble(
    WLE_rel = wle$WLE.rel %>% unique(), 
    EAP_rel = mod$EAP.rel
  )
}



###'#######################################################################'
###'
###' `gen_vec_beta()`
###' 
###' - Generate THREE different item difficulty vectors 
###'   based on the simulated true_theta 
###'
###' (1) `even`: evenly distributed values from min(theta) to max(theta)
###' (2) `uniform`: random samples from Unif(min(theta), max(theta))
###' (3) `normal`: random samples from N(mean = 0, sd = 1)
###'
###'

gen_vec_beta <- function(I = N_items, true_theta){
  
  # Extract summary statistics from the vector of true_theta
  min <- true_theta %>% min()
  max <- true_theta %>% max()
  
  # Simulate three different beta vectors
  beta_even <- c(0, seq(min, max, length = I - 1)) %>% sort()
  beta_unif <- runif(I, min = min, max = max) %>% sort()
  beta_norm <- rnorm(I, mean = 0, sd = 1) %>% sort()
  
  # Combine as a tibble (rowwise)
  df_temp <- tibble(c("beta_even", "beta_unif", "beta_norm"), 
                    list(beta_even, beta_unif, beta_norm)) %>%
    set_names(c("beta_kind", "vec_beta"))
  
  return(df_temp)
}




