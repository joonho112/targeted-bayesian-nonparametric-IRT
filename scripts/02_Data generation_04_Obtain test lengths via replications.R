
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



###'#######################################################################'
###'
###'  Generate a tibble crossing all person- and reliability-related 
###'  simulation factors 
###' 
###' 1) Number of persons (test takers)
###' 2) Data-generating functions for the latent trait (\theta) distribution
###' 3) true variance of \theta
###' 4) WLE reliability
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
df_theta_conds2
df_theta_conds2$true_theta[[30]]
df_theta_conds2$true_theta[[30]] %>% density() %>% plot() 



###'#######################################################################'
###'
###' Generate place holders for 100 beta replications
###'
###' - with the same \theta distributions 
###'   => the same `rep` number means the same `true_theta`
###'
###'

df_theta_conds3 <- df_theta_conds2 %>%
  mutate(sim_cond = 1:n()) %>%
  dplyr::select(sim_cond, everything()) %>%
  slice(rep(1:n(), each = 100)) %>%
  group_by(sim_cond) %>%
  mutate(rep = 1:n()) %>%
  ungroup()

df_theta_conds3 %>%
  dplyr::select(sim_cond, DGM, N_person, true_theta, rep)

df_theta_conds3$rep



###'#######################################################################'
###'
###' Add a vector for the test length (I, number of items) candidates
###'
###' -> result in 75 conditions by 100 replications by 70 test length choices
###'    = 525,000 rows 
###'
###'

### Test length candidates
vec_I <- c(1:70)

df_theta_I <- df_theta_conds3 %>%
  mutate(I = rep(list(vec_I), nrow(.))) %>%
  unnest(I)


# ### Work only with the first simulation condition
# df_theta_I_sub <- df_theta_I %>%
#   filter(sim_cond <= 5)
df_theta_I_sub <- df_theta_I  # 525,000 rows



###'#######################################################################'
###'
###' Simulate item difficulty (\beta) vectors
###'
###' => 525,000 * 3 beta DGM = 1,575,000 rows
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


### Save and load the resulting tibble
file_path <- file.path(data_dir2, 
                       "df_with_simulated_beta.rds")

saveRDS(df_theta_I_beta, file_path)

df_theta_I_beta <- read_rds(file = file_path)


### Check the result
df_theta_I_beta %>%
  select(sim_cond, DGM, N_person, WLE_rel, true_theta, rep, I, beta_kind, vec_beta)



###'#######################################################################
###'
###' Split data to accomodate large computations
###' 
###' - Total 1,575,000 rows
###'
###'

### Define subgroups based on simulation conditions + 100 replications
df_theta_I_beta %>%
  count(sim_cond, DGM, N_person, WLE_rel, true_var, rep)


###' Nest the data with the defined subgroup factors
df_nest <- df_theta_I_beta %>%
  group_by(sim_cond, DGM, N_person, WLE_rel, true_var, rep) %>%
  nest() %>%
  mutate(
    cond_name = paste(
      sprintf("%02d", sim_cond), 
      DGM, 
      "N", sprintf("%03d", N_person), 
      "WLErel", sprintf("%02d", WLE_rel*100), 
      "rep", sprintf("%03d", rep), 
      sep = "_"
    )
  ) %>%
  relocate(cond_name, .before = data)

df_nest
df_nest$cond_name
df_nest$data[[1]]



###'#######################################################################
###'
###' `process_df_nest()`
###' Define a wrapper function to process the nested data
###'
###'

process_df_nest <- function(sim_cond, DGM, N_person, true_var, 
                            WLE_rel, rep, cond_name, data){
  
  ### Extract data & prepare subset data
  df_sub <- data %>%
    select(I, true_theta, vec_beta) %>%
    rename(
      N_item = I, 
      theta = true_theta, 
      beta = vec_beta) %>%
    mutate(N_person = map_dbl(.x = theta, .f = length)) %>%
    relocate(N_person, .before = N_item)
  
  ### Simulate item response data
  list_sim_data <- df_sub %>%
    pmap(gen_rasch_logit_prob_y)
  
  df_sim <- df_sub %>%
    mutate(
      logit = map(.x = list_sim_data, .f = 1), 
      prob = map(.x = list_sim_data, .f = 2), 
      y = map(.x = list_sim_data, .f = 3)
    )
  
  ### Estimate test information
  df_sim2 <- df_sim %>%
    mutate(
      test_info_est = map_dbl(.x = logit, 
                              .f = get_test_info_est)
    )
  
  ### Fit the Rasch model
  tic()
  
  df_sim3 <- df_sim2 %>%
    mutate(
      mod_fit = map(.x = y, 
                    .f = safe_fit_rasch)
    )
  
  toc()
  
  ### Split result and error
  df_sim4 <- df_sim3 %>%
    mutate(
      error = map(.x = mod_fit, .f = "error"), 
      result = map(.x = mod_fit, .f = "result"), 
      error_lgl = map_lgl(.x = error, .f = not_null)
    ) %>%
    select(-mod_fit)
  
  ### Calculate reliabilities 
  df_sim5 <- df_sim4 %>%
    mutate(
      rel_est_temp = map(.x = result, 
                         .f = safe_get_reliabilities), 
      rel_est = map(.x = rel_est_temp, .f = "result"), 
      WLE_rel_est = map(.x = rel_est, .f = "WLE_rel") %>%
        map_dbl(.f = ~if_else(is_null(.), NA_real_, .)), 
      EAP_rel_est = map(.x = rel_est, .f = "EAP_rel") %>%
        map_dbl(.f = ~if_else(is_null(.), NA_real_, .))
    ) %>%
    select(-rel_est_temp, -rel_est)
  
  ### Calculate biases
  test_info_true <- data$test_info %>% unique()
  WLE_rel_true <- WLE_rel
  EAP_rel_true <- data$EAP_rel %>% unique()
  
  df_sim6 <- df_sim5 %>% 
    mutate(test_info_bias = test_info_est - test_info_true, 
           WLE_rel_bias = WLE_rel_est - WLE_rel_true, 
           EAP_rel_bias = EAP_rel_est - EAP_rel_true)
  
  ### Filter out results based on WLE estimates (biases)
  df_sim_final <- df_sim6 %>%
    mutate(beta_kind = data$beta_kind) %>%
    relocate(beta_kind, .after = N_item) %>%
    
    # group by beta_kind
    group_by(beta_kind) %>%
    
    # # reduce rows: drop errors 
    # filter(error_lgl == FALSE) %>%
    
    # # keep only by 25 percentile
    # filter(abs(WLE_rel_bias) <= abs(WLE_rel_bias) %>% quantile(prob = 0.25)) %>%
    
    # reduce columns: drop variables
    select(-error, -result)
  
  ### Save the resulting tibble as .rds file
  file_path <- file.path(data_dir2, 
                         paste0(cond_name, ".rds"))
  
  saveRDS(df_sim_final, file_path)
  
}

### Validation
safe_fit_rasch <- safely(fit_rasch) 
not_null <- negate(is_null)
safe_get_reliabilities <- safely(get_reliabilities)

df_nest[1:10,] %>%
  pmap(process_df_nest)



###'#######################################################################
###'
###' Loop over 7,500 rows
###'
###'

## Prepare parallel computing
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 20)


### Let's roll!
tic()

df_nest %>%
  future_pmap(.f = process_df_nest, 
              .options = furrr_options(seed = NULL,
                                       chunk_size = NULL,
                                       scheduling = 10),
              .progress = TRUE)

toc()



###'#######################################################################
###'
###' Preliminary analysis to collect beta solutions
###' from the resulting data frame
###'
###'

###' Data containing working directory
###' 7500 files
###' Total 142GB
data_dir3 <- file.path(data_dir2, 
                       "processed_df_nest")


### Investigate the resulting data frame example
temp_path <- file.path(data_dir3,
                       "01_Gaussian_N_020_WLErel_50_rep_001.rds")

df_temp <- read_rds(temp_path)


### Extract only scalars
df_scalar <- df_temp %>%
  select(N_item, beta_kind, 
         ends_with("_est"), ends_with("_bias"), 
         error_lgl)


###' Extract only the cases with the least WLE reliability biases
###' group_by beta_kind, drop all error cases
df_min <- df_temp %>%
  filter(error_lgl == FALSE) %>%
  group_by(beta_kind) %>%
  mutate(
    min_WLE_bias = min(abs(WLE_rel_bias))
  ) %>%
  ungroup() %>%
  filter(min_WLE_bias == abs(WLE_rel_bias))

df_min


### Mutate file location with the condition name
df_nest <- df_nest %>%
  mutate(
    rds_path = file.path(data_dir3, 
                         paste0(cond_name, ".rds"))
  ) %>%
  relocate(rds_path, .after = cond_name)



###'#######################################################################
###'
###' `extract_I_solutions()`
###' 
###' Extract minimum WLE reliability bias solution
###'
###'

extract_I_solutions <- function(rds_path){
  
  # Import the .rds file
  df_temp <- read_rds(rds_path)
  
  #' Extract only the cases with the least WLE reliability biases
  #' group_by beta_kind, drop all error cases
  df_min <- df_temp %>%
    filter(error_lgl == FALSE) %>%
    group_by(beta_kind) %>%
    mutate(
      min_WLE_bias = min(abs(WLE_rel_bias))
    ) %>%
    ungroup() %>%
    filter(min_WLE_bias == abs(WLE_rel_bias))
  
  # Delete unnecessary columns & return results
  df_min %>%
    select(-N_person)
}


###' Apply the function across the 7500 rows
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 20)

tic()

df_nest_min <- df_nest %>%
  ungroup() %>%
  select(-data) %>%
  # slice(1:100) %>%  # slice for testing
  mutate(
    df_min = 
      future_map(
        .x = .$rds_path, 
        .f = extract_I_solutions, 
        .options = furrr_options(seed = NULL,
                                 chunk_size = NULL,
                                 scheduling = 10),
        .progress = TRUE
      )
  ) %>%
  select(-rds_path)
  
toc()



###'#######################################################################
###'
###' Tidy up the resulting `df_nest_min` and save as a .rds file
###'
###'

### Check file size: 3.2 GB
object.size(df_nest_min) %>% format(units = "GB")


### Unnest the data frame
df_min_unnest <- df_nest_min %>%
  unnest(df_min)


### Drop unnecessary columns
names(df_min_unnest)

df_min_unnest$error_lgl %>% table()

df_min_unnest <- df_min_unnest %>%
  select(-error_lgl)


### Save as a .rds file 
object.size(df_min_unnest) %>% format(units = "GB")

temp_path <- file.path(data_dir2, 
                       "df_N_items_solutions.rds")

saveRDS(df_min_unnest, temp_path)



###'#######################################################################
###'
###' Obtain test lengths (the number of items) that lead to 
###' the preset level of WLE reliability
###' via 100 replications
###'
###'

### Subset only scalars; remove list columns
df_min_sum <- df_min_unnest %>%
  select(-theta, -beta, -logit, -prob, -y)


### Summarize based on beta_kind
df_min_solution <- df_min_sum %>%
  group_by(sim_cond, DGM, N_person, true_var, WLE_rel, beta_kind) %>%
  summarize(
    N_item_mean = mean(N_item), 
    WLE_rel_bias_mean = mean(min_WLE_bias), 
    WLE_rel_est_mean = mean(WLE_rel_est)
  )

temp_path <- file.path(work_dir, 
                       "tables", 
                       "df_N_item_solutions_Final.csv")

write_csv(df_min_solution, temp_path)


### Which beta_kind does lead to the least required number of items?
df_beta_compare <- df_min_solution %>%
  group_by(N_person, WLE_rel, beta_kind) %>%
  summarize(
    N_item_mean = mean(N_item_mean), 
    WLE_rel_bias_mean = mean(WLE_rel_bias_mean), 
    WLE_rel_est_mean = mean(WLE_rel_est_mean)
  )

temp_path <- file.path(work_dir, 
                       "tables", 
                       "df_N_item_solutions_Final_beta_comparison.csv")

write_csv(df_beta_compare, temp_path)


