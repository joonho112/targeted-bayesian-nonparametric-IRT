
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Results
###' 
###' Task: Get site-level estimates from posterior sample data
###'       
###' Data: Posterior sample data (.rds files)
###' 
###' Date: 
###' - 2022-06-23: initiated
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



###'#######################################################################
###'
###' Set things up to process posterior samples
###'
###'

### Set data containing directories
path_Gaussian <- file.path(data_dir2, 
                           "posterior_sample_Gaussian_NIMBLE_rep_001-020")

path_DPinform <- file.path(data_dir2, 
                           "posterior_sample_DPinform_NIMBLE_rep_001-005")

path_DPdiffuse <- file.path(data_dir2, 
                           "posterior_sample_DPdiffuse_NIMBLE_rep_001-005")

folder_path <- path_Gaussian # for test


### Define a function to prepare the posterior sample processing
prepare_postsamp_process <- function(df_pre, folder_path, 
                                     model_name = "Gaussian", 
                                     remove_unmatched = TRUE){
  
  df_temp <- list.files(folder_path) %>%
    tibble() %>%
    set_names("file_name") %>%
    mutate(cond_name = str_remove(file_name, ".rds")) %>%
    mutate(file_path = file.path(folder_path, file_name))
  
  df_merged <- df_pre %>%
    full_join_track(df_temp, by = c("cond_name")) %>%
    mutate(model = model_name) %>%
    select(model, everything())
  
  if (remove_unmatched == TRUE){
    
    df_merged %>%
      filter(!is.na(file_path))
    
  } else if (remove_unmatched == FALSE){   
    
    df_merged
  } 
}


### Apply to each postsamp folder
df_Gaussian <- prepare_postsamp_process(df_pre, path_Gaussian, "Gaussian", TRUE)
df_DPinform <- prepare_postsamp_process(df_pre, path_DPinform, "DP-inform", TRUE)
df_DPdiffuse <- prepare_postsamp_process(df_pre, path_DPdiffuse, "DP-diffuse", TRUE)


### Bind as one tibble
df_init <- bind_rows(df_Gaussian, df_DPinform, df_DPdiffuse) %>%
  mutate(model = factor(model, levels = c("Gaussian", "DP-inform", "DP-diffuse")))

object.size(df_init) %>% format("MB")



###'#######################################################################'
###'
###' Define the function `gen_site_estimates2()`
###' 
###' -> Generate posterior summary estimates for 
###' 
###'    the site-specific parameters (`theta_p`, `beta_i`, `hyperparam`)
###' 
###'    and merge true values and ML estimates
###'    
###' -> Why version `2`? 
###'    This function takes the `file_path` of the posterior samples,
###'    not the `posterior sample` themselves
###'    to save memory usage
###' 
###' -> Reference book from the multisite trial study
###'     
###' `n_j`: implied site size per site j
###' 
###' `se2_j`: sampling variance (SE^2) per site j
###'    
###' `tau_j_true`: true values of tau_j's
###' 
###' `tau_j_ML`: Observed or maximum likelihood (ML) estimates of tau_j's
###' 
###' `tau_j_PM`: Posterior mean (PM) estimates of tau_j's
###' 
###' `tau_j_PSD`: Posterior SD estimates of tau_j's
###' 
###' `tau_j_CB`: "Constrained Bayes (CB)" estimates of tau_j's
###'             by Ghosh (1992)
###' 
###' `tau_j_GR`: "Triple Goal (GR)" estimates of tau_j's
###'             using algorithm defined in Shen and Louis (1998)
###'
###' `R_bar`: Posterior means of ranks of tau_j's (1 = lowest)
###' 
###' `R_hat`: Integer ranks of tau_j's (`=rank(rbar)`)
###' 
###' 

gen_site_estimates2 <- function(file_path, theta, beta, model){
  
  # (1) Import the posterior sample .rds file
  postsamp <- read_rds(file_path) 
  
  df_postsamp <- postsamp %>% 
    as.data.frame() %>%
    tibble() %>%
    set_names(colnames(postsamp))
  
  
  # (2) Generate posterior summary estimates: theta (person parameters)
  df_theta <- df_postsamp %>%
    dplyr::select(contains("theta")) %>%
    as.matrix() %>%
    HETOP::triple_goal() %>%
    mutate(theta_true = theta) %>%
    relocate(theta_true, .before = "theta_pm") %>%
    rename(
      'person_id' = 'index',
      'theta_rbar' = 'rbar', 
      'theta_rhat' = 'rhat'
    ) %>%
    tibble()
  
  # (3) Generate posterior summary estimates: beta (item parameters)
  df_beta <- df_postsamp %>%
    dplyr::select(contains("beta")) %>%
    as.matrix() %>%
    HETOP::triple_goal() %>%
    mutate(theta_true = beta) %>%
    relocate(theta_true, .before = "theta_pm") %>%
    rename(
      'item_id' = 'index',
      'theta_rbar' = 'rbar', 
      'theta_rhat' = 'rhat'
    ) %>%
    tibble() %>%
    rename_with(~str_replace(.x, "theta", "beta"))
  
  # (4) Generate posterior summary estimates: hyperparameters
  if (model %in% c("Gaussian")){ 
    
    df_hyper <- df_postsamp %>%
      dplyr::select(-contains("theta")) %>%
      dplyr::select(-contains("beta")) %>%
      summarize(
        nu1_fix = mean(nu1),
        nu2_fix = mean(nu2), 
        s_tht_mean = mean(sqrt(s2_tht)), 
        s_tht_sd = sd(sqrt(s2_tht)), 
        s_tht_min = min(sqrt(s2_tht)), 
        s_tht_max = max(sqrt(s2_tht))
      )
      
  } else if (model %in% c("DP-inform", "DP-diffuse")){ 
    
    df_sub <- df_postsamp %>%
      dplyr::select(-contains("theta")) %>%
      dplyr::select(-contains("beta"))
    
    # alpha
    df_alpha <- df_sub %>%
      select(a, b, alpha) %>%
      summarize(
        a_fix = mean(a), 
        b_fix = mean(b),
        alpha_mean_drv = a_fix/b_fix, # theoretically derived alpha mean 
        alpha_mean_est = mean(alpha), 
        alpha_sd_drv = sqrt(a_fix)/b_fix, # theoretically derived alpha sd 
        alpha_sd_est = sd(alpha)
      )
    
    # mu_tilde
    df_mutilde <- df_sub %>%
      dplyr::select(contains("muTilde")) %>%
      map_dfr(summary, .id = "stat")
    
    # s2_tilde
    df_s2tilde <- df_sub %>%
      dplyr::select(contains("s2Tilde")) %>%
      map_dfr(summary, .id = "stat")
    
    # clustering behaviors
    df_zi <- df_sub %>%
      dplyr::select(contains("zi")) %>%
      pivot_longer(contains("zi"), names_to = "zi", values_to = "cluster_id") %>%
      arrange(zi, cluster_id) %>%
      count(cluster_id) %>%
      mutate(percent= 100*(n/sum(n)))
      
    # df_Nclust <- df_sub %>%
    #   dplyr::select(contains("zi")) %>%
    #   map(unique, .id = "zi") %>%
    #   map_int(length) 
    
    df_hyper <- list(df_alpha, df_zi, df_mutilde, df_s2tilde)
  } 
  
  # (5) Return results
  list(df_theta, df_beta, df_hyper)
}


# ### Safe version function
# safe_gen_site_estimates2 <- safely(gen_site_estimates2)
# not_null <- negate(is_null)


# ### Example application
# df_init_row <- df_init %>%
#   filter(model == "DP-inform") %>%
#   slice_sample(n = 1)


### Test the defined function
df_temp <- df_init %>%
  select(file_path, theta, beta, model) %>%
  slice_head(n = 10) %>%
  mutate(
    list_temp = pmap(.l = ., .f = gen_site_estimates2)
  ) %>%
  mutate(
    df_theta = map(.x = list_temp, .f = 1), 
    df_beta = map(.x = list_temp, .f = 2), 
    df_hyper = map(.x = list_temp, .f = 3)
  ) %>% 
  select(-list_temp, -theta, -beta)
  
df_merged <- df_init %>%
  right_join(df_temp, by = c("file_path", "model"))



###'#######################################################################
###'
###' Apply the `gen_site_estimates2()` function to all rows
###'
###'

### Prepare parallel computation: Set the number of workers
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 10)


### Processed tibble
df_init

df_init %>% 
  count(model)


### Let's roll!
tic()

df_est <- df_init %>%
  select(file_path, theta, beta, model) %>%
  # slice_head(n = 10) %>%
  mutate(
    list_temp = future_pmap(
      .l = ., 
      .f = gen_site_estimates2, 
      .options = furrr_options(seed = NULL,
                               chunk_size = 100,
                               scheduling = 1),
      .progress = TRUE
    )
  )

toc()

object.size(df_est) %>% format("MB")


### Post-processing
df_est2 <- df_est %>%
  mutate(
    df_theta = map(.x = list_temp, .f = 1), 
    df_beta = map(.x = list_temp, .f = 2), 
    df_hyper = map(.x = list_temp, .f = 3)
  ) %>% 
  select(-list_temp, -theta, -beta)

rm(df_est)


### Merge into the original data
df_site_est <- df_init %>%
  right_join(df_est2, by = c("file_path", "model"))

object.size(df_site_est) %>% format("MB")


### Save the resulting dataset
save_path <- file.path(data_dir2, "df_site_est_temp.rds")

write_rds(df_site_est, save_path)


