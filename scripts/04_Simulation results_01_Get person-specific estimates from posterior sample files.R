
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


### Extract the example data
y_mat <- df_pre$y[[1000]]

df_pre %>% 
  count(N_person)



###'#######################################################################
###'
###' Example file
###'
###'

###' Define the file path containing the posterior samples: An Example
temp_path <- file.path(data_dir2, 
                       "posterior_sample_DP-diffuse_20_50_100_200")

model_name <- c("DP-diffuse")


###' Pick one case
list_files <- list.files(temp_path)

file_name <- sample(list_files, 1)

file_path <- file.path(temp_path, file_name)

idx <- which(df_pre$cond_name == file_name %>% str_remove(".rds"))  

df_pre_row <- df_pre[idx, ]


###'#######################################################################'
###'
###' Define the function `gen_site_estimates2()`
###' 
###' -> Generate posterior summary estimates for 
###' 
###'    the site-specific parameters (`tau_j`)
###' 
###'    and merge true values and ML estimates
###'    
###' -> Why version `2`? 
###'    This function takes the `file_path` of the posterior samples,
###'    not the `posterior sample` themselves
###'    to save memory usage
###' 
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

gen_site_estimates2 <- function(df_pre_row, 
                                file_path){
  
  #' (1) Import the posterior sample file
  if (str_detect(file_path, ".fst")){ 
    
    postsamp <- read_fst(file_path)
    
  } else if (str_detect(file_path, ".rds")){
    
    postsamp <- readRDS(file_path)
    
  }
  
  postsamp_theta <- postsamp[[1]]
  postsamp_hyper <- postsamp[[2]]
  postsamp_beta <- postsamp[[3]]
  
  
  # (2-1) Generate posterior summary estimates: theta (person parameters)
  df_est_theta <- postsamp_theta %>%
    dplyr::select(contains("theta_")) %>%
    as.matrix() %>%
    HETOP::triple_goal() %>%
    mutate(theta_true = df_pre_row$theta[[1]]) %>%
    relocate(theta_true, .before = "theta_pm") %>%
    tibble()
  
  # (2-2) Generate posterior summary estimates: beta (item parameters)
  df_est_beta0 <- postsamp_beta %>%
    dplyr::select(contains("beta")) %>%
    as.matrix() %>%
    HETOP::triple_goal() %>%
    mutate(beta_true = df_pre_row$beta[[1]][-1])
  
  temp_names <- names(df_est_beta0) %>%
    str_replace("theta_", "beta_")
  
  df_est_beta <- df_est_beta0 %>% 
    set_names(temp_names) %>%
    tibble()
  
  
  # (3) Merge with original simulation data
  df_merged <- sim_data %>%
    bind_cols(df_est) %>%
    dplyr::select(site_id = index, 
                  n_j, se2_j, 
                  tau_j_true = tau_j, 
                  tau_j_ML = tau_j_hat, 
                  tau_j_PM = theta_pm, 
                  tau_j_PSD = theta_psd, 
                  tau_j_CB = theta_cb, 
                  tau_j_GR = theta_gr, 
                  R_bar = rbar, 
                  R_hat = rhat)
  
  return(df_merged)
}

safe_gen_site_estimates2 <- safely(gen_site_estimates2)
not_null <- negate(is_null)


# ### Testing: apply the gen_site_estimates() function
# df_temp <- df_sims_sub3 %>%
#   slice(1:5) %>%
#   mutate(site_est = map2(.x = sim_data, 
#                          .y = file_path, 
#                          .f = gen_site_estimates2))



df_est_theta_long <- df_est_theta %>%
  select(-rbar, -rhat, -theta_psd) %>%
  pivot_longer(theta_true:theta_gr, names_to = "type", values_to = "estimate") %>%
  mutate(type = factor(type))

ggplot(data = df_est_theta_long, aes(x = estimate)) +
  geom_density() +
  facet_wrap(~type)



###'#######################################################################'
###'
###' Define a function: `collect_postsamp()`
###' 
###' To collect posterior samples as list-columns
###' 
###' 

###' Define the file path containing the posterior samples: An Example
temp_path <- file.path(data_dir2, 
                       "posterior_sample_DP-diffuse_20_50_100_200")

model_name <- c("DP-diffuse")


### Define a function to create a tibble with the `postsamp` list-column
collect_postsamp <- function(temp_path, model_name){
  
  # Import all simulation results as a list
  setwd(temp_path)
  
  list_posterior <- list.files(temp_path) %>%
    # head() %>% 
    map(readRDS) 

  # Convert to a tibble dataframe
  df_posterior <- tibble(
    
    cond_name = list.files(temp_path) %>%
      head() %>%
      str_remove(".rds"), 
    
    model = model_name, 
  
    postsamp_theta = list_posterior %>% map(.f = 1), 
    postsamp_hyper = list_posterior %>% map(.f = 2), 
    postsamp_beta = list_posterior %>% map(.f = 3)
    
  )
  
  return(df_posterior)
}



###'#######################################################################
###'
###' Collect posterior samples
###'
###'

###' (1) DP-diffuse
temp_path <- file.path(data_dir2, 
                       "posterior_sample_DP-diffuse_20_50_100_200")

model_name <- c("DP-diffuse")

df_posterior1 <- collect_postsamp(temp_path, model_name = "DP-diffuse")













