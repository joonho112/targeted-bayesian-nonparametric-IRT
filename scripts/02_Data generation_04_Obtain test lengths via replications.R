
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Data generation 
###' 
###' Task: Obtain test lengths (the number of items) that lead to 
###'       the preset level of WLE reliability
###'       via 100 replications
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2021-12-12: created
###' - 2022-01-23: finalized
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
library(kableExtra)

library(psych)
library(sirt)
library(TAM)
library(LaplacesDemon)
library(sn)

library(future) 
library(furrr)
library(progressr)
library(tictoc)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)


# ### Prepare parallel computation: Set the number of workers
# parallelly::availableCores()
# parallelly::availableWorkers()
# plan(multisession, workers = 10)



###'#######################################################################'
###'
###' Generate simulation conditions 
###' 
###' - Varying the WLE reliability
###'
###'

### Target WLE reliabilities and fixed true variances
WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
true_var <- c(1, 1, 1, 1, 1)


### Generate simulation conditions
df_simconds <- gen_sim_conds(true_var, WLE_rel)

round(df_simconds, 2) %>%
  kbl(caption = "Simulation conditions #02") %>%
  kable_minimal(full_width = F, html_font = "Cambria")



###'#######################################################################'
###'
###'  Generate a tibble crossing all `varying but fixed` 
###'  person- and reliability-related simulation factors 
###' 
###' - Number of persons (test takers)
###' - Data-generating functions for the latent trait (\theta) distribution
###' - true variance of \theta
###' - WLE reliability
###' 
###'

###' [Full version] Define vectors of the varying but fixed simulation factors
vec_N_persons <- c(20, 50, 100, 200, 500)
vec_DGMs <- c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")
vec_WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
vec_true_var <- c(1, 1, 1, 1, 1)


###' [Reduced version] Define vectors of the varying but fixed simulation factors
vec_N_persons <- c(20, 50, 100, 200, 500)
vec_DGMs <- c("Gaussian", "ALD", "Mixed")
vec_WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
vec_true_var <- c(1, 1, 1, 1, 1)


### Generate a tibble crossing 5*6*5 simulation factors
df_simconds <- tibble_sim_conditions(vec_N_persons, 
                                     vec_DGMs, 
                                     vec_true_var, 
                                     vec_WLE_rel)

df_simconds



###'#######################################################################'
###'
###' Generate latent traits \theta's
###'
###' - Generate a column that includes the simulated theta's
###'
###' 

df_theta_conds2 <- df_simconds %>%
  mutate(true_theta = map2(.x = DGM, 
                           .y = N_person, 
                           .f = ~gen_vec_theta(.x, .y)))

### Validation
df_theta_conds2$true_theta
df_theta_conds2$true_theta[[30]] %>% density() %>% plot() 



###'#######################################################################'
###'
###' Generate place holders for 100 beta replications
###'
###' - with the same \theta distributions 
###'
###'

df_theta_conds3 <- df_theta_conds2 %>%
  mutate(sim_cond = 1:n()) %>%
  select(sim_cond, everything()) %>%
  slice(rep(1:n(), each = 100)) %>%
  group_by(sim_cond) %>%
  mutate(rep = 1:n()) %>%
  ungroup()



###'#######################################################################'
###'
###' Add a vector for the test length (I, number of items) candidates
###'
###' -> result in 75 conditions by 100 replications by 70 test length choices
###'    = 525,000 rows 
###'
###'

vec_I <- c(1:70)

df_theta_I <- df_theta_conds3 %>%
  mutate(I = rep(list(vec_I), nrow(.))) %>%
  unnest(I)


# ### Work only with the first simulation condition
# df_theta_I_sub <- df_theta_I %>%
#   filter(sim_cond <= 5)
# 
# df_theta_I_sub  # Only 35,000 cases

df_theta_I_sub <- df_theta_I



###'#######################################################################'
###'
###' Simulate item difficulty (\beta) vectors
###'
###'

### Prepare parallel computation: Set the number of workers
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 20)


### Let's roll!
tic()

df_theta_I_beta <- df_theta_I_sub %>%
  # slice(1:10) %>% # slices for testing
  mutate(
    
    ###' STEP 1. Simulate item difficulty (\beta) vectors
    beta = future_map2(.x = I, 
                       .y = true_theta, 
                       .f = ~gen_vec_beta(.x, .y), 
                       .options = furrr_options(seed = NULL, 
                                                chunk_size = NULL, 
                                                scheduling = 10), 
                       .progress = TRUE)
  ) %>%
  unnest(beta)

toc()


### Check the result
df_theta_I_beta

df_theta_I_beta$beta_kind

df_theta_I_beta$vec_beta[[100]]



###'#######################################################################'
###'
###' Simulate item response data
###'
###'

### Obtain standalone data list
tic()

list_sim_data <- df_theta_I_sub %>%
  select(N_person, I, true_theta, vec_beta) %>%
  rename(
    N_item = I, 
    theta = true_theta, 
    beta = vec_beta) %>%
  pmap(gen_rasch_logit_prob_y)

toc()


### Save as .rds file
setwd(data_dir2)
saveRDS(list_sim_data, file = "list_sim_data.rds")



###'#######################################################################'
###'
###' Estimate test information
###' 
###' 

### Extract only logit matrix
list_logit_mat <- map(.x = list_sim_data, .f = 1)


### Define a function to estimate test information
get_test_info_est <- function(logit_mat){
  
  # Item information vector
  logit_item <- apply(logit_mat, 2, mean)
  
  # Test information scalar (estimated)
  prob_item <- exp(logit_item)/(1 + exp(logit_item))
  item_info <- prob_item * (1 - prob_item)
  test_info_est <- sum(item_info)
  return(test_info_est)
}


### Obtain test information estimates
tic()

vec_test_info_est <- map_dbl(.x = list_logit_mat, 
                             .f = get_test_info_est)

toc()



###'#######################################################################'
###'
###' Fit the Rasch model
###'
###'

### Extract only the item response matrix
list_y_mat <- map(.x = list_sim_data, .f = 3)


### Define a function to safely fit Rasch model
safe_fit_rasch <- safely(fit_rasch)


### Fit the model across 105,000 rows
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 12)

tic()

list_mod_fit <- future_map(.x = list_y_mat, 
                           .f = safe_fit_rasch, 
                           .options = furrr_options(seed = NULL,
                                                    chunk_size = NULL,
                                                    scheduling = 5),
                           .progress = TRUE)

toc()


### Save as .rds file
setwd(data_dir2)
saveRDS(list_mod_fit, file = "list_mod_fit.rds")



###'#######################################################################'
###'
###' Calculate reliabilities
###'
###'

### Extract only the result 
list_mod_result <- map(.x = list_mod_fit, .f = "result")


### Define a safe version of function to calculate reliabilities
safe_get_reliabilities <- safely(get_reliabilities)


# ### Estimate reliabilities (future_map version)
# tic()
# 
# list_rel_est <- future_map(.x = list_mod_result, 
#                            .f = safe_get_reliabilities, 
#                            .options = furrr_options(seed = NULL,
#                                                     chunk_size = NULL,
#                                                     scheduling = 10),
#                            .progress = TRUE) %>%
#   map(.f = "result")
# 
# toc()


### Estimate reliabilities (single-core version)
tic()

list_rel_est <- map(.x = list_mod_result, 
                    .f = safe_get_reliabilities) %>%
  map(.f = "result")

toc()


### Extract WLE and EAP reliability estimates
vec_WLE_est <- map(.x = list_rel_est, .f = "WLE_rel") %>%
  map_dbl(.f = ~if_else(is_null(.), NA_real_, .))

vec_EAP_est <- map(.x = list_rel_est, .f = "EAP_rel") %>%
  map_dbl(.f = ~if_else(is_null(.), NA_real_, .))



###'#######################################################################'
###'
###' Attach to the original dataframe 
###' and calculate biases between estimates and true fixed values
###' 
###' (1) test_info_est - test_info
###' (2) WLE_rel_est - WLE_rel
###' (3) EAP_rel_est - EAP_rel  
###' 
###' 

###' Attach to the original dataframe
df_theta_I_sub2 <- df_theta_I_sub %>%
  mutate(
    sim_data = list_sim_data, 
    test_info_est = vec_test_info_est, 
    WLE_rel_est = vec_WLE_est, 
    EAP_rel_est = vec_EAP_est
  )


###' Calculate biases
df_theta_I_sub3 <- df_theta_I_sub2 %>%
  mutate(bias_test_info = test_info_est - test_info, 
         bias_WLE_rel = WLE_rel_est - WLE_rel, 
         bias_EAP_rel = EAP_rel_est - EAP_rel)



###'#######################################################################'
###'
###' Filter out the best results based on WLE estimates (biases)
###'
###'

### Filter out cases with minimum biase WLE reliability
df_solution <- df_theta_I_sub3 %>%
  group_by(sim_cond, rep, beta_kind) %>%
  filter(abs(bias_WLE_rel) == min(abs(bias_WLE_rel), na.rm = TRUE))


### Select variables to present
vec_vars <- c("sim_cond", "rep",
              "beta_kind", "DGM", "N_person", "I", 
              "test_info", "test_info_est", 
              "WLE_rel", "WLE_rel_est", 
              "EAP_rel", "EAP_rel_est", 
              "bias_test_info", "bias_WLE_rel", "bias_EAP_rel")

df_solution2 <- df_solution %>%
  select(all_of(vec_vars))


### Average over 100 replications
df_solution3 <- df_solution2 %>%
  group_by(sim_cond, beta_kind, DGM, N_person, 
           test_info, WLE_rel, EAP_rel) %>%
  summarize(mean_I = mean(I, na.rm = TRUE), 
            
            mean_WLE_rel_est = mean(WLE_rel_est, na.rm = TRUE), 
            mean_EAP_rel_est = mean(EAP_rel_est, na.rm = TRUE), 
            mean_test_info_est = mean(test_info_est, na.rm = TRUE),
            
            mean_bias_WLE_rel = mean(bias_WLE_rel, na.rm = TRUE), 
            mean_bias_EAP_rel = mean(bias_EAP_rel, na.rm = TRUE), 
            mean_bias_test_info = mean(bias_test_info, na.rm = TRUE))
            

### Save the results
setwd(work_dir)
write_csv(df_solution3, file = "tables/test length solutions.csv")
