###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Data generation 
###' 
###' Task: Define a preliminary function 
###'       to simulate item response data according to
###'       the level of reliabilities
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2021-10-25
###' - 2021-11-14 
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
library(kableExtra)
library(psych)
library(sirt)
library(TAM)
library(LaplacesDemon)
library(sn)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################'
###'
###' Generate simulation conditions
###' 
###' (1) Varying the separation coefficients
###'
###'

### True variance and SD
true_var <- c(1, 1, 1, 1, 1)
true_SD <- sqrt(true_var)


###' Signal-to-noise ratio
###' = Separation coefficient
###' = true SD / error SD (RMSE)
sep_coef <- c(1, 2, 3, 4, 5)


### Error variance and RMSE
err_SD <- true_SD/sep_coef    
err_var <- err_SD^2

### Observed variance = true variance + error variance
obs_var <- true_var + err_var


###' (Person separation) reliability = true variance / observed variance
###' EAP reliability is larger than WLE reliability
WLE_rel <- true_var/obs_var  # or, 1 - (err_var/obs_var)
EAP_rel <- 1 - (err_var/(obs_var + err_var))


### Strata = (4*sep_coef + 1)/3
strata <- (4*sep_coef + 1)/3


### Combine as a tibble
tab_cond <- tibble(
  true_SD, true_var, err_SD, err_var, obs_var, 
  sep_coef, strata, WLE_rel, EAP_rel
)

round(tab_cond, 2) %>%
  kbl(caption = "Simulation conditions #01") %>%
  kable_minimal(full_width = F, html_font = "Cambria")



###'#######################################################################'
###'
###' Generate simulation conditions
###' 
###' (2) Varying the WLE reliability
###'
###'

# Target WLE reliabilities and fixed true variances
WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
true_var <- c(1, 1, 1, 1, 1)

# Define a function to generate a table with simulation conditions
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
  
  # EAP reliability
  EAP_rel <- 1 - (err_var/(obs_var + err_var))
  
  # Combine as a tibble
  tab_cond <- tibble(
    true_SD, true_var, err_SD, err_var, obs_var, 
    sep_coef, strata, WLE_rel, EAP_rel
  )
  
  return(tab_cond)
}

# New simulation conditions
df_simconds <- gen_sim_conds(true_var, WLE_rel)
round(df_simconds, 2) %>%
  kbl(caption = "Simulation conditions #02") %>%
  kable_minimal(full_width = F, html_font = "Cambria")



###'#######################################################################'
###'
###' Generating theta's - true latent trait distributions
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


### Simulate theta's according to the defined DGMs
N <- 10000

Gaussian <- gen_true_theta(N = N, true_dist = "Gaussian", tau = 0, var = 1)

T <- gen_true_theta(N = N, true_dist = "T", tau = 0, var = 1, nu = 5)

Skew <- gen_true_theta(N = N, true_dist = "Skew", tau = 0, var = 1, slant = 5)

ALD <-gen_true_theta(N = N, true_dist = "ALD", tau = 0, var = 1, rho = 0.1)

Bimodal <- gen_true_theta(N = N, true_dist = "Mixture", tau = 0, var = 1,
                          delta = 4, eps = 0.3, ups = 1)

Mixed <- gen_true_theta(N = N, true_dist = "Mixture", tau = 0, var = 1,
                        delta = 5, eps = 0.3, ups = 2)

### Collect the simulated data
vec_levels <- c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")

df_theta <- tibble(Gaussian, T, Skew, ALD, Bimodal, Mixed) %>%
  pivot_longer(Gaussian:Mixed, names_to = "DGM", values_to = "value") %>%
  mutate(DGM = factor(DGM, levels = vec_levels)) %>%
  arrange(DGM)


### Plot 
p <- ggplot(data = df_theta, aes(x = value)) + 
  geom_density() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  facet_wrap(DGM ~ .) + 
  xlim(-6, 6) + 
  theme_bw() + 
  labs(title = "Data-generating models for true person ability distributions", 
       subtitle = "with mean 0 and variance 1", 
       x = NULL, y = "Density", 
       caption = "Range truncated at [-6, 6]")

file_path <- file.path(work_dir, 
                       "figures", 
                       "True person ability distributions.pdf")

ggsave(file_path, p, width = 8, height = 5)



###'#######################################################################'
###'
###' Generating a vector of SEM
###' 
###' (1) A simple hierarchical approach 
###' 
###'     - assuming the same SEMs across the theta range
###'
###'     - Doesn't work well
###'     
###'

### Generate true theta's from a Gaussian distribution
N <- 100
true_theta <- gen_true_theta(N = N, true_dist = "Gaussian", tau = 0, var = 1)

tibble(true_theta) %>%
  round(2) %>%
  arrange(true_theta) %>%
  pull(true_theta)


### Extract SEMs from the pre-defined simulation conditions
err_SD <- df_simconds %>% pull(err_SD)


### Extract the first SEM and generate a vector with a length of N
vec_err_SD1 <- rep(err_SD[2], N)


### Level-1 sampling model
obs_theta <- rnorm(n = N, mean = true_theta, sd = vec_err_SD1)


### Present the resulting true theta, observed theta, and SEM table
tibble(true_theta, obs_theta, vec_err_SD1) %>%
  set_names("true", "observed", "SEM") %>%
  arrange(true) %>%
  print(n = 50)


# Calculate the target number of items with Dr. Linacre's approximation 
# and some variations
N_item <- prophecy_target_N_item(C = 20, R_C = 0.88, R_T = 0.6) %>% round()


# Define a function to generate binary responses
gen_rasch_prob_y <- function(N_person, N_item, theta, beta){
  
  # Create two empty matrices saving probabilities and responses 
  y_mat <- prob_mat <- eta_mat <- matrix(0, nrow = N_person, ncol = N_item) 
  
  # Generate responses based on Rasch model
  for(i in 1:N_person) { 
    for(j in 1:N_item) { 
      
      eta <- theta[i] - beta[j]
      eta_mat[i, j] <- eta
      
      prob <- exp(eta)/(1 + exp(eta))
      prob_mat[i, j] <- prob
      
      y_mat[i,j] <- rbinom(1, 1, prob) 
    } 
  }
  
  # Return three matrices as a list
  dimnames(y_mat) <- dimnames(prob_mat) <- dimnames(eta_mat) <-
    list(
      paste0("id_", str_pad(seq(1:N_person), 2, pad = "0")), 
      paste0("item_", str_pad(seq(1:N_item), 2, pad = "0"))
    )
  list(eta_mat, prob_mat, y_mat)
}


### Example simulation
N <- 100
I <- N_item

theta <- obs_theta
beta <- c(0, seq(-2, 2, length = I - 1)) %>% sort()

list_prob_y <- gen_rasch_prob_y(N, I, theta, beta)

list_prob_y


y_mat <- list_prob_y[[3]]




# Fit a simple Rasch model (MML estimation) 
mod <- TAM::tam.mml(resp = list_prob_y[[2]])


# Extract item parameters 
mod$item  # item difficulties


# WLE estimation
wle <- TAM::tam.wle(mod)

# Standard errors
se <- TAM::tam.se(mod)

se$xsi


# EAP, WLE reliabilities
mod$EAP.rel

wle$WLE.rel %>% unique()



###'#######################################################################'
###'
###' Generating a vector of SEM
###' 
###' (2) Using the test information function
###'    
###'     - assuming different SEMs across the theta range
###'
###'

### Calculating the test information functions
err_SD <- df_simconds %>% 
  pull(err_SD)

(SEM <- err_SD)  # standard error of measurement

(test_info <- 1/(SEM^2))  # test information function


### Target test information level = 9
test_info[1]



# Set the number of items
# Sum of the 20 item information values = 9
# on average, 9/20 = 0.45
I <- 9
beta <- runif(I, min = -2, max = 2) %>% sort()


# Calculate the probabilities of correct answer
list_prob_y <- gen_rasch_prob_y(N, I, true_theta, beta)
eta_mat <- list_prob_y[[1]]

# Marginalize out person dimension
eta_item <- apply(eta_mat, 2, mean)

# Get a item information vector
prob_item <- exp(eta_item)/(1 + exp(eta_item))
item_info <- prob_item * (1- prob_item)
sum(item_info)

mod <- TAM::tam.mml(resp = list_prob_y[[3]])

# Extract item parameters 
mod$item  # item difficulties


# WLE estimation
wle <- TAM::tam.wle(mod)

# Standard errors
se <- TAM::tam.se(mod)

se$xsi


# EAP, WLE reliabilities
mod$EAP.rel

wle$WLE.rel %>% unique()
