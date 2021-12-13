
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Data generation 
###' 
###' Task: Simulate data according to the level of reliabilities
###'       - Using test information function
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2021-11-15 
###' - 2021-11-28 finalized
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
data_dir2 <- c("~/Documents/Data-files/targeted-bayesian-nonparametric-IRT-large-files")
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
###' Generate latent traits theta's
###' 
###' - All possible SIX DGMs:
###'   Gaussian, T, Skew, ALD, Bimodal, Mixed
###'   
###' - Varying number of persons    
###'
###'

### Varying number of persons
vec_N_persons <- c(20, 50, 100, 200, 1000)


### Varying data-generating functions
vec_DGMs <- c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")


### A table for all combinations 
df_theta_conds <- list(
  DGM = vec_DGMs, 
  N_person = vec_N_persons
) %>% cross_df()


### Attach df_simconds to each row
df_theta_conds <- df_theta_conds %>%
  mutate(sim_cond = rep(list(df_simconds), nrow(df_theta_conds))) 


### Generate a column that includes the simulated theta's
df_theta_conds2 <- df_theta_conds %>%
  mutate(true_theta = map2(.x = DGM, 
                           .y = N_person, 
                           .f = ~gen_vec_theta(.x, .y)))

# df_theta_conds2$true_theta[[28]] %>% density() %>% plot() # valiation


### Unnest simulation conditions
df_theta <- df_theta_conds2 %>%
  unnest(sim_cond)


### Save the simulation conditions as a table
df_tab_simconds <- df_theta %>% dplyr::select(-true_theta)

file_path <- file.path(work_dir, 
                       "tables", 
                       "simulation_conditions.csv")

write_csv(df_tab_simconds, file_path)



###'#######################################################################'
###'
###' Generating a vector of SEM
###' 
###' - Using the test information function
###'   assuming different SEMs across the theta range
###'
###'

### Calculate the target level of test information
df_theta <- df_theta %>%
  mutate(test_info = 1/err_var) %>%
  relocate(test_info, .before = WLE_rel)


### Loop over 150 conditions and save the results
list_grid_search <- list()  # empty list to store grid search results
list_solution <- list()     # empty list to store solutions from grid search

for (i in seq(nrow(df_theta))){
  
  #' A vector of N_items for grid search
  #' grid search option - Number of items = {1, ..., 99, 100}
  df_I <-  c(1:100)%>% 
    tibble() %>%
    set_names("I")
  
  # Slice a row, attach df_I, and expand it
  df_row <- df_theta[i, ] %>%
    mutate(df_I = list(df_I)) %>%
    unnest(df_I)
  
  #' Generate item difficulty vectors 
  #' and expand the dataframe by three beta vectors
  df_row2 <- df_row %>%
    mutate(
      beta = map2(.x = I, .y = true_theta, 
                  .f = ~gen_vec_beta(.x, .y))
    ) %>% 
    unnest(beta)
  
  #' (1) Simulate data
  list_sim_data <- df_row2 %>%
    select(N_person, I, true_theta, vec_beta) %>%
    rename(
      N_item = I, 
      theta = true_theta, 
      beta = vec_beta) %>%
    pmap(gen_rasch_logit_prob_y)
  
  #' (2) Estimate test information
  list_logit_mat <- map(.x = list_sim_data, .f = 1)
  
  get_test_info_est <- function(logit_mat){
    
    # Item information vector
    logit_item <- apply(logit_mat, 2, mean)
    
    # Test information scalar (estimated)
    prob_item <- exp(logit_item)/(1 + exp(logit_item))
    item_info <- prob_item * (1 - prob_item)
    test_info_est <- sum(item_info)
    return(test_info_est)
  }
  
  vec_test_info_est <- map_dbl(.x = list_logit_mat, 
                               .f = get_test_info_est)
  
  #' (3) Fit the Rasch model
  list_y_mat <- map(.x = list_sim_data, .f = 3)
  
  safe_fit_rasch <- safely(fit_rasch)
  
  lit_mod_fit <- map(.x = list_y_mat, 
                     .f = safe_fit_rasch)
  
  #' (4) Calculate reliabilities
  list_mod_result <- map(.x = list_mod_fit, .f = "result")
  
  safe_get_reliabilities <- safely(get_reliabilities)
  
  list_rel_est <- map(.x = list_mod_result, .f = safe_get_reliabilities) %>%
    map(.f = "result")
  
  vec_WLE_est <- map(.x = list_rel_est, .f = "WLE_rel") %>%
    map_dbl(.f = ~if_else(is_null(.), NA_real_, .))
  
  vec_EAP_est <- map(.x = list_rel_est, .f = "EAP_rel") %>%
    map_dbl(.f = ~if_else(is_null(.), NA_real_, .))
  
  
  #' Attach to df_row2
  df_row3 <- df_row2 %>%
    mutate(
      sim_data = list_sim_data, 
      test_info_est = vec_test_info_est, 
      WLE_rel_est = vec_WLE_est, 
      EAP_rel_est = vec_EAP_est
      )
  
  #' Calculate the biases between
  #' (1) test_info_est - test_info
  #' (2) WLE_rel_est - WLE_rel
  #' (3) EAP_rel_est - EAP_rel  
  df_grid_search <- df_row3 %>%
    mutate(bias_test_info = test_info_est - test_info, 
           bias_WLE_rel = WLE_rel_est - WLE_rel, 
           bias_EAP_rel = EAP_rel_est - EAP_rel)
  
  list_grid_search[[i]] <- df_grid_search
  
  # Filter out the best result based on WLE estimates
  df_solution <- df_grid_search %>% 
    group_by(beta_kind) %>%
    filter(abs(bias_WLE_rel) == min(abs(bias_WLE_rel), na.rm = TRUE))
  
  list_solution[[i]] <- df_solution
}



###'#######################################################################'
###'
###' Save the resulting SEM vector and simulated data
###'
###'

### Embed grid search results and solutions to df_theta
df_result <- df_theta %>%
  mutate(
    grid_search = list_grid_search, 
    solution = list_solution
  )

df_result


### Obtain soultions as vectors
subset_vars <- function(df){
  
  vec_vars <- c("beta_kind", "DGM", "N_person", "I", 
                "test_info", "test_info_est", 
                "WLE_rel", "WLE_rel_est", 
                "EAP_rel", "EAP_rel_est")
  
  df %>%
    ungroup() %>%
    select(all_of(vec_vars)) %>%
    arrange(beta_kind)
    
}

df_result <- df_result %>%
  mutate(
    solution_sub = map(.x = solution, 
                       .f = subset_vars)
  )


### Save the results
setwd(data_dir2)
saveRDS(df_result, file = "simulation_results.rds")
