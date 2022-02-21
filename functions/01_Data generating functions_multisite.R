
###'######################################################################
###'
###' Project: Multisite Variation Project
###'
###' Category: Define functions
###' 
###' Task: Define functions for data generation
###'
###' - Updated data-generating functions with three-step approach
###' - 1) Generate the prior distribution G (with zero mean and unit variance)   
###' - 2) Generate the sampling (or standard) errors
###' - 3) Generate observed site-specific effects 
###'      
###' Date: 
###' 2022-01-05  created
###' 2022-01-09  updated 
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###'

###'######################################################################
###'
###' Call necessary libraries
###'
###'

library(tidyverse)
library(LaplacesDemon)
library(sn)



###'#######################################################################'
###'
###' [STEP 1] `gen_priorG()`
###'    
###' - A function to generate prior distribution G's 
###'   with zero mean and unit variance
###'     
###' - [updated] Rescale this G with cross-site variation
###' 
###' (1) Gaussian
###' 
###' (2) Student T distribution with df = `nu`
###' 
###' (3) Skewed normal distribution 
###'     with marginal mean and variance of SkewN(0, 1^2)
###' 
###' (3) Asymmetric Laplace distribution with skewness parameter `p`
###' 
###' (4) A mixture of two Gaussian distributions
###'    - `delta`: distance between two means
###'    - `eps`: proportion of the small component
###'    - `ups`: ratio between two variances
###'   
###'   * Bimodal: delta = 4; eps = 0.3; ups = 1
###'   * Mixed: delta = 5; eps = 0.3; ups = 2
###'
###'

gen_priorG <- function(true_dist = NULL, 
                       J = 50, sigma_tau = 0.25, 
                       nu = NULL, slant = NULL, rho = NULL, 
                       delta = NULL, eps = NULL, ups = NULL){

  # Set zero mean and unit variance for all Gs
  tau <- 0; sigma_tau_unit <- 1L; var = 1L 
  
  # Generate true site-specific effects tau_j
  
  if (true_dist == "Gaussian"){ 
    
    tau_j <- rnorm(J, mean = tau, sd = sigma_tau_unit)
    
  } else if (true_dist == "T" & !is.null(nu)){   
    
    tau_j <- rt(J, nu)*sqrt((nu - 2)/nu)
    
  } else if (true_dist == "Skew" & !is.null(slant)){
    
    # Generate location and scale parameters of skewed normal
    # to achieve E(theta) = 0 and Var(theta) = 1
    delta <- slant/sqrt(1 + slant^2)
    scale <- sqrt(var/(1-(2*delta^2/pi)))
    location <- tau - scale*sqrt(2/pi)*delta
    
    tau_j <- sn::rsn(n = J, xi = location, omega = scale, 
                     alpha = slant)[1:J]
    
  } else if (true_dist == "ALD" & !is.null(rho)){
    
    # Generate location, scale, and skewness parameters 
    # of ALD distribution to achieve E(tau_j) = 0 and Var(tau_j) = 1
    scale <- sqrt((2*rho^2*var)/(1 + rho^4))
    location <- tau - ((scale*(1/rho - rho))/sqrt(2))
    tau_j <- LaplacesDemon::ralaplace(J, location, scale, rho)

  } else if (true_dist == "Mixture" & !(is.null(delta)|is.null(eps)|is.null(ups))){
    
    # Define a normalizing factor `a`
    a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
    
    # Simulate a mixture of two normals with mean 0 and variance 1  
    ind <- runif(J) < (1 - eps)
    tau_j <- ind*rnorm(J, -eps*delta/a, sqrt(1/a^2)) + 
      (1 - ind)*rnorm(J, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  }
  
  ###' Rescale the resulting tau_j
  ###' unit SD * cross-site impact SD
  tau_j_scaled <- tau_j*sigma_tau
  
  return(tau_j_scaled)
}


# ### Check out the function 
# list_tau_j <- rerun(.n = 1000, 
#                     gen_priorG(J = 25, true_dist = "Gaussian", sigma_tau = 0.05))
# 
# list_tau_j <- rerun(.n = 1000, 
#                     gen_priorG(J = 25, true_dist = "T", sigma_tau = 0.05, 
#                                nu = 5))
# 
# list_tau_j <- rerun(.n = 1000, 
#                     gen_priorG(J = 25, true_dist = "Skew", sigma_tau = 0.05, 
#                                slant = 5))
# 
# list_tau_j <- rerun(.n = 1000, 
#                     gen_priorG(J = 160, true_dist = "ALD", sigma_tau = 0.25, 
#                                rho = 0.1))
# 
# list_tau_j <- rerun(.n = 1000, 
#                     gen_priorG(J = 160, true_dist = "Mixture", sigma_tau = 0.25, 
#                                delta = 4, eps = 0.3, ups = 1))
# 
# list_tau_j <- rerun(.n = 1000, 
#                     gen_priorG(J = 160, true_dist = "Mixture", sigma_tau = 0.25, 
#                                delta = 5, eps = 0.3, ups = 2))
# 
# map_dbl(.x = list_tau_j, .f = mean) %>% mean()
# map_dbl(.x = list_tau_j, .f = sd) %>% mean()



###'#######################################################################'
###'
###' [STEP 2] `gen_nj_se2j_vec_gamma()`
###'
###' - A function to generate a vector of 
###'   simulated within-site sampling errors (sampling variances)
###'   based on the simulated site sizes
###'   
###' - [updated] generate the site size vector with gamma distribution
###'   continuous gamma distribution with rounding + censoring at nj_min = 5
###' 
###' (1) `J`: number of sites {25, 50, 75, 100, 300}
###' (2) `nj_mean`: mean site sizes {10, 20, 40, 80, 160}
###' (3) `cv`: coefficient of variation. CV = sd/mean {0, 0.25, 0.5, 0.75}
###' (4) `nj_min`: lower bound of the possible site sizes
###' 
###' (5) `p`: proportion of units treated
###' (6) `R2`: proportion explained by covariates
###' 
###' 

gen_nj_se2j_vec_gamma <- function(J = 25, 
                                  nj_mean = 10, 
                                  cv = 0, 
                                  nj_min = 5, 
                                  p = 0.50, 
                                  R2 = 0.00){
  
  # (1) Simulate site sizes based on a Gamma distribution
  if (cv == 0){ 
    
    nj_vec <- rep(nj_mean, J)
    
  } else if (cv != 0){   
    
    # Simulate continuous Gamma distribution
    a <- 1/cv^2  
    b <- a/nj_mean
    nj_raw_gamma <- rgamma(n = J, shape = a, rate = b)
    
    # Censor at the preset lower bound nj_min = 5
    # Round to nearest integer
    nj_vec <- if_else(nj_raw_gamma <= nj_min, nj_min, nj_raw_gamma) %>%
      round(0)
  }
  
  # (2) Calculate SE^2 from the simulated site size vector
  varY <- 1L  # in effect size units
  kappa <- varY * (1/p + 1/(1-p)) * (1 - R2)
  se2_vec <- kappa * (1/nj_vec)  # scaled with kappa
  
  # Return both: site sizes & standard errors (scale: variance, not SD)
  df_se2 <- data.frame(nj_vec, se2_vec)
  names(df_se2) <- c("n_j", "se2_j")
  return(df_se2)
}


# ### Check out the function 
# df_se2 <- gen_nj_se2j_vec_gamma(J = 50, nj_mean = 20, cv = 0.75, 
#                                 p = 0.50, R2 = 0.00)



###'#######################################################################'
###'
###' [STEP 3] `gen_tau_j_hat()`
###' 
###' - A function to generate observed site-specific estimates: tau_j_hat
###' 
###' (1) tau_j's: Prior distribution G of true site-specific effects
###' (2) n_j's: simulated site sizes
###' (3) se2_j's: within-site sampling "variances", not SDs
###' (4) tau_j_hat's: estimated/observed site-specific effects
###'
###'

# ### Obtain tau_j
# tau_j <- gen_priorG(J = 50, true_dist = "Mixture", sigma_tau = 0.25, 
#                     delta = 5, eps = 0.3, ups = 2)
# G_dist(tau_j)
# 
# 
# ### Obtain n_j and SE_j
# df_se2 <- gen_nj_se2j_vec_gamma(J = 50, nj_mean = 20, cv = 0.75,
#                                 p = 0.50, R2 = 0.00)

gen_tau_j_hat <- function(tau_j = NULL, df_se2 = NULL){
  
  df_temp <- df_se2 %>% 
    # Shuffle the order of rows (n_j, se2_j)
    slice_sample(prop = 1) %>% 
    
    # Combine true site-specific effects
    mutate(tau_j = tau_j) %>%  
    
    # Generate a vector of observed site-specific effects
    mutate(tau_j_hat = rnorm(n = length(tau_j), 
                             mean = tau_j, 
                             sd = sqrt(se2_j))) %>%
    
    # Reorder rows by site sizes n_j
    arrange(n_j)
  
  return(df_temp)
}


# ### Check out the function
# gen_tau_j_hat(tau_j, df_se2)

