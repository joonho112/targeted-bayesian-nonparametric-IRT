
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Implementation
###' 
###' Task: Estimate Gaussian, DP-diffuse, DP-inform models using NIMBLE
###'       
###' Data: Posterior sample data (.rds files)
###' 
###' Date: 
###' - 2022-06-25
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

library(HETOP)
library(nimble)

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
    N_person == 100, 
    WLE_rel == 0.9, 
    rep == 5
  )

Y <- df_example$y[[1]]



###'#######################################################################
###'
###' Define handy functions
###'
###'

### A function to complie NIMBLE function
compile_and_run_NIMBLE <- function(nimble_code, constants, data, inits, 
                                   niter = 4000, nburnin = 2000){
  
  # Create and compile model
  model <- nimbleModel(nimble_code, constants, data, inits)
  
  cModel <- compileNimble(model)
  
  
  ### Build an MCMC to fit the model
  conf <- configureMCMC(model, monitors = monitors)
  
  model_MCMC <- buildMCMC(conf)
  
  cModel_MCMC <- compileNimble(model_MCMC, project = model)
  
  
  ### Let's roll! - Run an MCMC and obtain posterior samples
  postsamp <- runMCMC(cModel_MCMC, niter = niter, nburnin = nburnin)
  
  return(postsamp)
}



### A function for quickly checking site-specific estimation
check_site_specific_results <- function(postsamp = postsamp_Gaussian1, 
                                        param = "theta", 
                                        df_example = df_example){
  
  #' Generate site-specific estimates
  #' Merge with true values
  df_est <- postsamp %>%
    data.frame() %>%
    dplyr::select(contains(param)) %>%
    as.matrix() %>%
    HETOP::triple_goal() %>%
    tibble() %>%
    rename_with(.fn = ~str_replace(., "theta", param), 
                .cols = everything()) %>%
    mutate(!!paste0(param, "_true") := df_example[[param]][[1]]) %>%
    relocate(!!paste0(param, "_true"), .after = "index")
  
  
  # A quick plot for comparing estimates
  df_long <- df_est %>%
    dplyr::select(index, contains(param), -contains("_psd")) %>% 
    pivot_longer(cols = contains(param), 
                 names_to = "estimator",
                 values_to = "estimate")
  
  p <- ggplot(data = df_long, 
              aes(x = estimate, color = estimator, group = estimator)) +
    geom_density() 

  # Return the resulting dataframe and plot
  list(df_est, p)
}


### DPM model checker
DPM_model_check <- function(postsamp){
  
  # (1) alpha ~ gamma(a, b) parameters
  params_alpha <- postsamp %>% 
    as.data.frame() %>%
    select(all_of(c("a", "b"))) %>%
    as.matrix() %>%
    samplesSummary()
  
  # (2) alpha posterior distribution
  p1 <- postsamp %>% 
    as.data.frame() %>%
    select("alpha") %>%
    ggplot(aes(x = alpha)) + 
    geom_density() +
    theme_minimal() + 
    labs(title = "posterior distribution of alpha (concentration parameter)")
  
  # (3) muTilde
  p2 <- postsamp %>%
    as.data.frame() %>%
    dplyr::select(contains("muTilde")) %>% 
    map(.f = mean) %>%
    unlist() %>%
    tibble() %>%
    set_names("muTilde") %>%
    ggplot(aes(x = muTilde)) + 
    geom_density() +
    theme_minimal() + 
    labs(title = "posterior distribution of N muTilde's")
  
  # (4) s2Tilde
  p3 <- postsamp %>%
    as.data.frame() %>%
    dplyr::select(contains("s2Tilde")) %>% 
    map(.f = mean) %>%
    unlist() %>%
    tibble() %>%
    set_names("s2Tilde") %>%
    ggplot(aes(x = s2Tilde)) + 
    geom_density() +
    theme_minimal() + 
    labs(title = "posterior distribution of N s2Tilde's")
  
  # (5) zi N vectors
  p_zi <- postsamp %>%
    as.data.frame() %>%
    select(contains("zi")) %>% 
    pivot_longer(everything(), names_to = "zi", values_to = "est") %>%
    mutate(zi = factor(zi)) %>%
    ggplot(aes(x = est, group = zi)) + 
    geom_density(size = 0.1) + 
    theme(legend.position = "none") + theme_bw() + 
    labs(title = "The CRP-distributed vector zi (N overlaid lines)", 
         x = "cluster ID")
  
  
  p_hyper <- plot_grid(p1, p2, p3, nrow = 3,
                       labels = "AUTO", label_size = 12, 
                       align = "v")  
  
  # Return results
  list(params_alpha, p_hyper, p_zi)
  
}



###'#######################################################################
###'
###' Estimate the Bayesian Rasch model
###' 
###' (1) Gaussian prior 
###'    - #1. Fixed sigma_theta (theta[p] ~ dnorm(0, 1))
###'
###'

### Define NIMBLE model code
code_Gaussian1 <- nimbleCode({ 
  
  for(i in 1:I) {
    for(p in 1:N) {
      y[p,i] ~ dbern(pi[p,i])
      logit(pi[p,i]) <-  theta[p] - beta[i]
    }
    beta[i] ~ dnorm(0, var = 10)
  }  
  
  for(p in 1:N) {
    theta[p] ~ dnorm(0, 1)
  }
  
})


### Set data, constants, monitors, and inits
data <- list(y = Y)

constants <- list(I = ncol(Y), N = nrow(Y))

monitors <- c("beta", "theta")

inits <- list(
  beta = rnorm(constants$I, 0, 1), 
  theta = rnorm(constants$N, 0, 1)
)


### Create, compile, build, and run MCMC model
tic()

postsamp_Gaussian1 <- compile_and_run_NIMBLE(nimble_code = code_Gaussian1, 
                                             constants, data, inits, 
                                             niter = 4000, nburnin = 2000)

toc()


### Check and tidy up the results
postsamp_Gaussian1 %>% 
  dim()

(summary_Gaussian1 <- samplesSummary(postsamp_Gaussian1))


### Check site-specific estimates
(list_est_theta <- check_site_specific_results(postsamp = postsamp_Gaussian1,
                                              param = "theta", 
                                              df_example = df_example))

(list_est_beta <- check_site_specific_results(postsamp = postsamp_Gaussian1,
                                             param = "beta", 
                                             df_example = df_example))



###'#######################################################################
###'
###' Estimate the Bayesian Rasch model
###' 
###' (2) Gaussian prior 
###'    - #2. Hyperprior on sigma_theta
###'
###'

### Define NIMBLE model code
code_Gaussian2 <- nimbleCode({ 
  
  for(i in 1:I) {
    for(p in 1:N) {
      y[p,i] ~ dbern(pi[p,i])
      logit(pi[p,i]) <-  theta[p] - beta[i]
    }
  }  
  
  for(i in 1:I) {
    beta[i] ~ dnorm(0,  10)
  } 
  
  for(p in 1:N) {
    theta[p] ~ dnorm(0, var = s2_tht)
  }  
  
  # mu_bt ~ dnorm(0, var = 3)
  # s2_bt ~ dinvgamma(2.01, 1.01)
  
  s2_tht ~ dinvgamma(2.01, 1.01)
  
})


### Set data, constants, monitors, and inits
data <- list(y = Y)

constants <- list(I = ncol(Y), N = nrow(Y))

monitors <- c("beta", "theta", 
              # "mu_bt", "s2_bt", 
              "s2_tht")

inits <- list(
  beta = rnorm(constants$I, 0, 1), 
  theta = rnorm(constants$N, 0, 1), 
  # mu_bt = 0, s2_bt = 1, 
  s2_tht = 1
)


### Create, compile, build, and run MCMC model
tic()

postsamp_Gaussian2 <- compile_and_run_NIMBLE(nimble_code = code_Gaussian2, 
                                             constants, data, inits, 
                                             niter = 4000, nburnin = 2000)

toc()


### Check and tidy up the results
postsamp_Gaussian2 %>% 
  dim()

(summary_Gaussian2 <- samplesSummary(postsamp_Gaussian2))


### Check site-specific estimates
(list_est_theta <- check_site_specific_results(postsamp = postsamp_Gaussian2,
                                               param = "theta", 
                                               df_example = df_example))

(list_est_beta <- check_site_specific_results(postsamp = postsamp_Gaussian2,
                                              param = "beta", 
                                              df_example = df_example))



###'#######################################################################
###'
###' Estimate the Bayesian Rasch model
###' 
###' (3) Dirichlet Process Mixture model
###'    - #1. DP-basic model
###'    -     The model from Paganin et al. (2021)
###'
###'

### Define NIMBLE model code
code_DPbasic <- nimbleCode({ 
  
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


### Set data, constants, monitors, and inits
data <- list(y = Y)

constants <- list(I = ncol(Y), N = nrow(Y), M = 50)

monitors <- c("beta", "theta", 
              "zi", "alpha", 
              "a", "b", 
              "muTilde", "s2Tilde")

inits <- list(
  beta = rnorm(constants$I, 0, 1), 
  alpha = 1, a = 1, b = 3, 
  nu1 = 2.01, nu2 = 1.01, s2_mu = 2
)

scores <- apply(data$y, 1, sum)
std_scores <- (scores - mean(scores))/sd(scores)
inits$theta <- std_scores
inits$zi <- kmeans(std_scores, 3)$cluster


### Create, compile, build, and run MCMC model
tic()

postsamp_DPbasic <- compile_and_run_NIMBLE(nimble_code = code_DPbasic, 
                                           constants, data, inits, 
                                           niter = 4000, nburnin = 2000)

toc()


### Check and tidy up the results
postsamp_DPbasic %>% 
  dim()

(summary_DPbasic <- samplesSummary(postsamp_DPbasic))


### Check site-specific estimates
(list_est_theta <- check_site_specific_results(postsamp = postsamp_DPbasic,
                                               param = "theta", 
                                               df_example = df_example))

(list_est_beta <- check_site_specific_results(postsamp = postsamp_DPbasic,
                                              param = "beta", 
                                              df_example = df_example))


### DPM model check
DPbasic_check <- DPM_model_check(postsamp = postsamp_DPbasic)

DPbasic_check[1]
DPbasic_check[2]
DPbasic_check[3]



###'#######################################################################
###'
###' Estimate the Bayesian Rasch model
###' 
###' (4) Dirichlet Process Mixture model
###'    - #2. DP-diffuse model
###'
###'

### Define NIMBLE model code
code_DPdiffuse <- nimbleCode({ 
  
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
  
  # for(m in 1:M) {
  #   muTilde[m] ~ dnorm(0, var = 200)
  #   s2Tilde[m] ~ dinvgamma(1, 1)
  # }
  
})


### Set data, constants, monitors, and inits
data <- list(y = Y)

constants <- list(I = ncol(Y), N = nrow(Y), M = 500)

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


### Create, compile, build, and run MCMC model
tic()

postsamp_DPdiffuse <- compile_and_run_NIMBLE(nimble_code = code_DPdiffuse, 
                                           constants, data, inits, 
                                           niter = 4000, nburnin = 2000)

toc()


### Check and tidy up the results
postsamp_DPdiffuse %>% 
  dim()

(summary_DPdiffuse <- samplesSummary(postsamp_DPdiffuse))


### Check site-specific estimates
(list_est_theta <- check_site_specific_results(postsamp = postsamp_DPdiffuse,
                                               param = "theta", 
                                               df_example = df_example))

(list_est_beta <- check_site_specific_results(postsamp = postsamp_DPdiffuse,
                                              param = "beta", 
                                              df_example = df_example))


### DPM model check
DPdiffuse_check <- DPM_model_check(postsamp = postsamp_DPdiffuse)

DPdiffuse_check[1]
DPdiffuse_check[2]
DPdiffuse_check[3]



###'#######################################################################
###'
###' Estimate the Bayesian Rasch model
###' 
###' (5) Dirichlet Process Mixture model
###'    - #3. DP-inform model
###'
###'

### Define NIMBLE model code
code_DPinform <- nimbleCode({ 
  
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


### Set data, constants, monitors, and inits
data <- list(y = Y)

constants <- list(I = ncol(Y), N = nrow(Y), M = 250)

monitors <- c("beta", "theta", 
              "zi", "alpha", 
              "a", "b", 
              "muTilde", "s2Tilde")

a <- 9.32; b <- 0.88
alpha_mean <- nrow(Y)*0.1
alpha_var <- a/b^2  

inits <- list(
  beta = rnorm(constants$I, 0, 1), 
  alpha = alpha_mean, 
  a = a, b = b, 
  nu1 = 2.01, nu2 = 1.01, s2_mu = 2
)

scores <- apply(data$y, 1, sum)
std_scores <- (scores - mean(scores))/sd(scores)
inits$theta <- std_scores
inits$zi <- kmeans(std_scores, 5)$cluster


### Create, compile, build, and run MCMC model
tic()

postsamp_DPinform <- compile_and_run_NIMBLE(nimble_code = code_DPinform, 
                                             constants, data, inits, 
                                             niter = 4000, nburnin = 2000)

toc()


### Check and tidy up the results
postsamp_DPinform %>% 
  dim()

(summary_DPinform <- samplesSummary(postsamp_DPinform))


### Check site-specific estimates
(list_est_theta <- check_site_specific_results(postsamp = postsamp_DPinform,
                                               param = "theta", 
                                               df_example = df_example))

(list_est_beta <- check_site_specific_results(postsamp = postsamp_DPinform,
                                              param = "beta", 
                                              df_example = df_example))


### DPM model check
DPinform_check <- DPM_model_check(postsamp = postsamp_DPinform)

DPinform_check[1]
DPinform_check[2]
DPinform_check[3]



###'#######################################################################
###'
###' Compare three different Gaussian and DPM model specifications
###'
###' (1) Person parameters - theta
###'
###'

### Collect site-specific estimates: theta
list_theta_Gaussian1 <- check_site_specific_results(
  postsamp = postsamp_Gaussian1,
  param = "theta", 
  df_example = df_example
)

list_theta_Gaussian2 <- check_site_specific_results(
  postsamp = postsamp_Gaussian2,
  param = "theta", 
  df_example = df_example
)

list_theta_DPbasic <- check_site_specific_results(
  postsamp = postsamp_DPbasic,
  param = "theta", 
  df_example = df_example
)

list_theta_DPdiffuse <- check_site_specific_results(
  postsamp = postsamp_DPdiffuse,
  param = "theta", 
  df_example = df_example
)  

list_theta_DPinform <- check_site_specific_results(
  postsamp = postsamp_DPinform,
  param = "theta", 
  df_example = df_example
)  

# plot_grid(list_theta_DPbasic[[2]], 
#           list_theta_DPdiffuse[[2]], 
#           list_theta_DPinform[[2]], 
#           nrow = 3,
#           labels = "AUTO", label_size = 12, 
#           align = "v")  


### Tidy up estimated thetas
df_theta_Gauss1 <- list_theta_Gaussian1[[1]] %>%
  mutate(model = "Gaussian1") %>%
  select(model, everything())

df_theta_Gauss2 <- list_theta_Gaussian2[[1]] %>%
  mutate(model = "Gaussian2") %>%
  select(model, everything())

df_theta_basic <- list_theta_DPbasic[[1]] %>%
  mutate(model = "DPM-basic") %>%
  select(model, everything())

df_theta_diffuse <- list_theta_DPdiffuse[[1]] %>%
  mutate(model = "DPM-diffuse") %>%
  select(model, everything())

df_theta_inform <- list_theta_DPinform[[1]] %>%
  mutate(model = "DPM-inform") %>%
  select(model, everything())

lev_model <- c("Gaussian1", "Gaussian2", 
               "DPM-basic", "DPM-diffuse", "DPM-inform")

### Make a long dataframe
df_theta_wide <- bind_rows(
  df_theta_Gauss1, 
  df_theta_Gauss2, 
  df_theta_basic, 
  df_theta_diffuse,
  df_theta_inform) %>%
  mutate(model = factor(model, levels = lev_model))

df_theta_long <- df_theta_wide %>%
  select(-theta_psd, -rbar, -rhat) %>%
  pivot_longer(contains("theta"), names_to = "estimator", values_to = "estimate") %>%
  mutate(estimator = str_remove(estimator, "theta_"), 
         estimator = str_to_upper(estimator), 
         estimator = factor(estimator, 
                            levels = c("TRUE", "PM", "CB", "GR")))
  

### Plot
ggplot(df_theta_long, aes(x = estimate, group = model, color = model)) + 
  geom_density() +
  facet_grid(rows = vars(estimator))


ggplot(df_theta_long, aes(x = estimate, group = estimator, color = estimator)) + 
  geom_density() +
  facet_grid(rows = vars(model))



