
###'######################################################################
###'
###' Category: Simulation implementation
###' 
###' Task: Estimating Dirichlet Process Mixture (DPM) Models
###'       
###'      - Select a prior for precision parameter alpha
###'        
###'      - Calculate a, b parameters for alpha ~ Gamma(a, b)
###'      
###'      - Apply my `Novel DORO method` using a grid search
###'      
###'      => Define functions for my `Novel DORO method`
###'      
###' Data: Simulated data
###' 
###' Date: 
###' - 2020-03-19: initiated
###' - 2021-09-04: updated
###' - 2021-12-05: update for defining functions
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###' 

###'######################################################################
###'
###' Basic settings
###'
###'

### Start with a clean slate
gc(); rm(list=ls())   


### Set working directory 
#work_dir <- c("~/Documents/Bayes-deconvolution")
work_dir <- c("~/targeted-bayesian-nonparametric-IRT") # for Windows
setwd(work_dir)


### Set a data directory
data_dir <- file.path(work_dir, "datasets")


### Call libraries
library(tidyverse)
library(cowplot)
library(bspmma)
library(rstan)
library(gmp)
library(broom)
library(metR)


### Call libraries for parallel computation
library(future) 
library(furrr)
library(progressr)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)


### Preset ggplot themes
theme_preset <- theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())



###'#######################################################################'
###'
###'  `get_induced_K_prior()`
###'  
###'  a function to generate the alpha induced prior of K
###'  `K = the assumed number of latent clusters`
###' 
###' - n: Number of random draws for the gamma distribution
###' - a: shape parameter of the gamma distribution
###' - b: rate parameter of the gamma distribution  
###' - J: Number of sites
###' 
###' 

get_induced_K_prior <- function(n, a, b, J){
  
  # Generate the assumed alpha vector
  alpha <- rgamma(n = n, shape = a, rate = b)
  
  # Define an empty vector to store pobability mass
  temp_vec <- vector()
  
  for (k in seq(J)){
    
    #' Encode the constant part using the
    #' "unsigned" Stirling number of the first kind
    sign <- (-1)^(J - k)
    constant <- sign*(b^a*as.numeric(gmp::Stirling1(n = J, k = k)))/gamma(a)
    
    # Encode the numerical integration part
    integrand <- function(alpha, J, k, a, b){
      
      log_comp1 <- (k + a - 1)*log(alpha)
      log_comp2 <- -b*alpha
      log_comp3 <- lgamma(alpha)
      log_comp4 <- -lgamma(alpha + J)
      
      ifelse(alpha == 0, 0, exp(log_comp1 + log_comp2 + log_comp3 + log_comp4))
    }
    
    integral <- integrate(f = integrand, lower = 0, upper = Inf, 
                          J = J, k = k, a = a, b = b)$value
    
    # Store in the empty vector
    temp_vec[k] <- constant*integral
  }
  
  return(temp_vec)
}


### Test the function
n <- 2000
a <- 1.78
b <- 2.24
J <- 35
induced_K_prior_ex <- get_induced_K_prior(n, a, b, J)
plot(density(induced_K_prior_ex))



###'#######################################################################'
###'
###' `get_df_K_priors()`
###' a function to generate the alpha induced prior of K 
###' `K = the assumed number of latent clusters`
###' along with other reference distributions
###' 
###' - n: Number of random draws for the gamma distribution
###' - a: shape parameter of the gamma distribution
###' - b: rate parameter of the gamma distribution  
###' - J: Number of sites
###' - dof: degrees of freedom of the chi-squared distribution
###' 
###'  1. `Uniform`: Uniform distribution of K
###'  2. `Chi-squared`: Chi-squared distribution of K with degrees of freedom
###'  3. `alpha-induced`: alpha(n, a, b)-induced prior of K
###' 
###'

get_df_K_priors <- function(n, a, b, J, dof){
  
  # Generate uniform, chi-squared, alpha-induced K distributions
  df_temp <- list(
    K = seq(J),  # K index
    Uniform = rep(1/J, J), # Uniform K dist. 
    Chi_squared = dchisq(seq(J), df = dof, ncp = 0), # Chi-squared K dist.
    alpha_induced = get_induced_K_prior(n, a, b, J) # alpha-induced K dist.
  ) %>%
    bind_cols()
  
  # Return a list along with the input parameters
  list(
    df_K_priors = df_temp, 
    n = n, a = a, b = b, J = J, dof = dof
    )
}

safe_get_df_K_priors <- safely(get_df_K_priors)

### Test the function 
n <- 2000
a <- 1.78
b <- 2.24
J <- 35

list_K_priors <- get_df_K_priors(n = 2000, a = 1.78, b = 1.50, J = 35, dof = 3.5)



###'#######################################################################'
###'
###' `plot_K_priors()`
###' a function to plot the alpha induced prior of K 
###' `K = the assumed number of latent clusters`
###' along with other reference distributions
###' 
###' - Take the `list_K_priors` as a input, and plot the following:
###' 
###'  1. `Uniform`: Uniform distribution of K
###'  2. `Chi-squared`: Chi-squared distribution of K with degrees of freedom
###'  3. `alpha-induced`: alpha(n, a, b)-induced prior of K
###' 
###'

plot_K_priors <- function(list_K_priors){
  
  # Prepare a dataframe to plot
  vec_dist <- c("Uniform", "Chi_squared", "alpha_induced")
  
  df_plot <- list_K_priors[["df_K_priors"]] %>%
    pivot_longer(Uniform:alpha_induced, 
                 names_to = "K_distribution", 
                 values_to = "value") %>%
    mutate(K_distribution = factor(K_distribution, 
                                   levels = vec_dist)) %>%
    arrange(K_distribution)
  
  # Prepare labels
  lab_chi2 <- paste0("Chi2 distribution with df = ", list_K_priors$dof)
  lab_alpha <- paste0("alpha-induced distribution with Gamma(", 
                      list_K_priors$a, ", ", list_K_priors$b, ")")
  
  labels_temp <- labs(
    title = "The probability mass function for the expected number of clusters (K)", 
    subtitle = paste0(lab_chi2, " vs. ", lab_alpha), 
    x = paste0("K: The expected number of clusters (J = ", list_K_priors$J, ")"), 
    y = "Probability"
    ) 
  
  # plot
  ggplot(data = df_plot, aes(x = K, 
                             y = value, 
                             group = K_distribution, 
                             color = K_distribution, 
                             linetype = K_distribution)) + 
    geom_line(size = 1) + 
    scale_x_continuous(breaks = seq(from = 0, to = J, by = 2)) + 
    labels_temp + theme_preset
}


### Test the function 
plot_K_priors(list_K_priors)



###'######################################################################
###'
###' `get_KLD_info()`
###' 
###' 
###' a function to calculate the Kullback-Leibler divergence measure
###' 
###' between (1) The prior of K and (2) the alpha-induced prior for K
###'
###' - pi_x: the “true” distribution of data, observations, 
###'         or theoretical distribution
###'         = Chi-squared(dof) -> encoded belief about the prior distribution of K
###' 
###' - pi_y: a theory, model, or approximation of pi_x
###'         = Gamma(a, b) -> precision parameter alpha
###'
###'

get_KLD_info <- function(pi_x, pi_y){
  
  # Get a logical index if pi_y's not NaN or Inf
  idx <- (!is.infinite(pi_y) & !is.nan(pi_y))
  
  # Subset only valid cases
  pi_y_sub <- pi_y[idx]
  pi_x_sub <- pi_x[idx]
  
  ### Return KLD measure
  KLD <- sum(pi_x_sub*(log(pi_x_sub) - log(pi_y_sub)))
  return(KLD)
}

### Test the function
df_K_priors <- list_K_priors[["df_K_priors"]]

safe_get_KLD_info <- safely(get_KLD_info) # safe version

safe_get_KLD_info(pi_x = df_K_priors$Chi_squared, 
                  pi_y = df_K_priors$alpha_induced) # test the function



###'#######################################################################'
###'
###' `ab_grid_search()`
###'  
###'  a function to search for the shape (`a`) and rate (`b`) parameters 
###'  for gamma distribution that makes the following two distributions
###'  as close as possible:
###'  
###'  1. `Chi-squared`: The prior distribution of K
###'     - encoded by the chi-squared distribution with dof = K
###'  
###'  2. `alpha-induced`: The alpha-induced distribution of K
###'     - the precision parameter alpha ~ Gamma(a, b) 
###'     
###' - `J`: Number of sites
###' - `dof`: degrees of freedom of the chi-squared distribution
###' 
###' Feed `J` and `dof`, and then get `a` and `b`
###'  
###'  

### Prepare parallel computation: Set the number of workers
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 20)


ab_grid_search <- function(J, dof, 
                           a_seq = c(from = 0, to = 20, by = 1),
                           b_seq = c(from = 0, to = 20, by = 1)){
  
  ###' Prepare wrapped functions
  safe_get_df_K_priors <- safely(get_df_K_priors)
  safe_get_KLD_info <- safely(get_KLD_info)
  safe_ab_grid_search <- safely(ab_grid_search)
  not_null <- negate(is_null)
  
  ###' (1) Set up an input grid
  ###'     - both a, b ranges from 0 to 10
  # n_breaks <- 501
  # x <- seq(0, 10, length = n_breaks)[-1]
  
  grid_a <- seq(a_seq[1], a_seq[2], by = a_seq[3])[-1]
  grid_b <- seq(b_seq[1], b_seq[2], by = b_seq[3])[-1]
  
  df_grid <- list(
    J = J, dof = dof, a = grid_a, b = grid_b
  ) %>% 
    cross_df() %>%
    mutate(n = 2000) %>%
    dplyr::select(n, everything())
  
  
  ###' (2) Calculate df_K_priors
  df_temp <- df_grid %>%
    # slice(1:10) %>%
    mutate(
      # df_K_priors list
      df_K_priors = 
        future_pmap(.l = ., 
                    .f = safe_get_df_K_priors, 
                    .options = furrr_options(seed = NULL),
                    .progress = TRUE), 
      # Extract result and error
      result = map(.x = df_K_priors, .f = "result"), 
      df_K = map(.x = result, .f = "df_K_priors"), 
      
      error_list = map(.x = df_K_priors, .f = "error"), 
      error_K = map_lgl(.x = error_list, .f = not_null)
    ) %>%
    dplyr::select(n, J, dof, a, b, result, df_K, error_K) %>%
    filter(error_K == FALSE) # sort out errors
  
  ###' (3) Calculate KLD measure
  df_KLD <- df_temp %>%
    mutate(
      KLD_list = map(.x = df_K, 
                     .f = ~safe_get_KLD_info(pi_x = .x$Chi_squared, 
                                             pi_y = .x$alpha_induced)), 
      KLD = map_dbl(.x = KLD_list, .f = "result"), 
      error_KLD = map(.x = KLD_list, .f = "error"), 
      error_KLD = map_lgl(.x = error_KLD, .f = not_null)
    ) %>%
    dplyr::select(-KLD_list)
  
  return(df_KLD)
}


### Test the function
tic()

df_ab_grid <- ab_grid_search(J = 150, dof = 50)

toc()

idx <- which.min(abs(df_ab_grid$KLD))

ab_solution <- df_ab_grid[idx, ]

plot_K_priors(ab_solution$result[[1]])



###'######################################################################
###'
###' `plot_ab_solution()`
###' 
###' Get and plot ab solution with KLD measure
###' 
###' 

plot_ab_solution <- function(df_ab_grid){
  
  ###' (1) Get a, b solution
  ###'     that minimized the KLD measure
  idx <- which.min(abs(df_ab_grid$KLD))
  
  ab_solution <- df_ab_grid[idx, ]
  
  p_ab_solution <- plot_K_priors(ab_solution$result[[1]]) # ab solution plot
  
  vec_ab <- c(ab_solution$a, ab_solution$b)
  
  J <- ab_solution$J
  dof <- ab_solution$dof
  
  ###' (2) Generate a dataframe to plot
  ###'     restrict a, b range as 3 times of each a, b solution
  df_plot <- df_ab_grid %>%
    filter(a < 3*ab_solution$a, b < 3*ab_solution$b)
    
  ###' (3) Define labels
  label_temp <-   labs(
    title = 
      paste0("The Kullback-Leibler divergence measure (J = ", J, ", df = ", dof, ")"), 
    subtitle = 
      "Assuming an informative prior on the number of clusters (K)", 
    x = "Shape parameter (a) of the Gamma prior", 
    y = "Rate parameter (b) of the Gamma prior", 
    caption = 
      paste0("Solution: ", 
             "a = ", ab_solution$a, ", ", 
             "b = ", ab_solution$b, ", ", 
             "KLD = ", round(ab_solution$KLD, 3))
    )
   
  ###' (4) Generate a contour plot
  p_raster <- ggplot(data = df_plot, aes(x = a, y = b, z = KLD)) +
    geom_raster(aes(fill = KLD)) +
    scale_fill_viridis_c(direction = 1) +  
    geom_contour(color = "white") +
    geom_text_contour(colour = "white") + 
    stat_subset(aes(subset = abs(KLD) < abs(min(KLD))*3), color = "red") + 
    theme_bw() + label_temp
  
  ### Return resulting plots
  list(ab_solution, p_ab_solution, p_raster)
}


### Test the function
list_solution <- plot_ab_solution(df_KLD)
list_solution[[2]]
list_solution[[3]]

