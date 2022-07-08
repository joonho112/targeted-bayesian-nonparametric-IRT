
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Data generation 
###' 
###' Task: Final data generation 
###'       
###' Data: Simulated data
###' 
###' Date: 
###' - 2022-07-07
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
                      # "Documents",
                      "targeted-bayesian-nonparametric-IRT") 

data_dir <- file.path(work_dir, "datasets")

data_dir2 <- file.path(path.expand("~"), 
                       # "Documents",
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

save_path <- file.path(work_dir, 
                       "tables", 
                       "simulation_conditions_summary.csv")

write_csv(df_simconds, save_path)




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
View(df_simconds)



###'#######################################################################'
###'
###' Expand the simulation conditions: for 100 replications
###' 
###'

df_simconds_by_reps <- df_simconds %>%
  mutate(sim_cond = 1:n()) %>%
  dplyr::select(sim_cond, everything()) %>%
  slice(rep(1:n(), each = 100)) %>%
  group_by(sim_cond) %>%
  mutate(rep = 1:n()) %>%
  ungroup()

df_simconds_by_reps  

View(df_simconds_by_reps)



###'#######################################################################'
###'
###' Generate latent traits \theta's
###'
###' - Generate a column that includes the simulated theta's
###' - I made a mistake at this phase in the previous analysis turn
###'
###' 

df_sim_theta <- df_simconds_by_reps %>%
  mutate(true_theta = map2(.x = DGM, 
                           .y = N_person, 
                           .f = ~gen_vec_theta(.x, .y)))

### Validation
View(df_sim_theta)

test_row <- df_sim_theta %>%
  filter(
    DGM %in% c("Mixed"),
    N_person %in% c(500), 
    WLE_rel %in% c(0.5), 
    rep == 11
) 

test_row$true_theta[[1]] %>% 
  density() %>%
  plot()



###'#######################################################################'
###'
###' Simulate item difficulty data 
###' 
###' via Monte Carlo approximation
###' 
###' (1) Add a vector for the test length (I, number of items) candidates
###'
###' -> result in 
###'    75 conditions by 
###'    100 replications by 
###'    70 test length choices
###'  
###'  = 525,000 rows 
###'
###'

### Test length candidates
vec_I <- c(1:70)

df_sim_theta_I <- df_sim_theta %>%
  mutate(I = rep(list(vec_I), nrow(.))) %>%
  unnest(I)

df_sim_theta_I

View(df_sim_theta_I)



###'#######################################################################'
###' 
###' Simulate item difficulty data 
###' 
###' via Monte Carlo approximation
###' 
###' (2) Simulate item difficulty (\beta) vectors
###'
###' => 525,000 * 3 beta DGM = 1,575,000 rows
###' 
###' Chop into `3*5*5 = 75` files for memory savings
###' 
###'

### Prepare chopped looping for memory savings
vec_DGM <- c("Gaussian", "ALD", "Mixed")
vec_N_person <- c(20, 50, 100, 200, 500)
vec_WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
vec_rep <- seq(100)


### Construct a for loop for memory savings
tic()

for (i in seq_along(vec_DGM)){
  for (j in seq_along(vec_N_person)){
    for (k in seq_along(vec_WLE_rel)){
      
      # Print progress
      cat(paste0("DGM: ", vec_DGM[i]), 
          paste0("|| N_person:", vec_N_person[j]), 
          paste0("|| WLE_rel: ", vec_WLE_rel[k]),
          "\n")
      
      # Prepare parallel computation: Set the number of workers
      # parallelly::availableCores()
      # parallelly::availableWorkers()
      plan(multisession, workers = 8)
      
      # Subset the conditions
      df_sim_sub <- df_sim_theta_I %>%
        filter(DGM %in% vec_DGM[i]) %>%
        filter(N_person %in% vec_N_person[j]) %>%
        filter(WLE_rel %in% vec_WLE_rel[k])
      
      # Let's roll!
      tic()
      
      df_temp <- df_sim_sub %>%
        # slice_head(n = 70000) %>% # slices for testing
        mutate(
          
          ###' STEP 1. Simulate item difficulty (\beta) vectors
          beta = future_map2(.x = I, 
                             .y = true_theta, 
                             .f = ~gen_vec_beta(.x, .y), 
                             .options = furrr_options(seed = NULL, 
                                                      chunk_size = 2000, 
                                                      scheduling = 5), 
                             .progress = TRUE)
        ) %>%
        unnest(beta)
      
      toc()
      
      # Save to the disk
      file_name <- paste0("df", 
                          "_DGM-", vec_DGM[i], 
                          "_N-", vec_N_person[j], 
                          "_WLE_rel-", vec_WLE_rel[k]*100, ".rds")
      
      save_path <- file.path(data_dir2, 
                             "temp_folder", 
                             file_name)
      
      write_rds(df_temp, save_path)
    }
  }
}

toc()



###'#######################################################################'
###' 
###' Simulate item difficulty data 
###' 
###' via Monte Carlo approximation
###' 
###' (3) Simulate `item response data` based on theta and beta 
###'
###' - Loop over `3*5*5 = 75` files for memory savings
###' 
###' - Work only on the `beta_norm` for now to save time
###' 
###'

### Prepare chopped looping for memory savings
vec_DGM <- c("Gaussian", "ALD", "Mixed")
vec_N_person <- c(20, 50, 100, 200, 500)
vec_WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# list_rep <- split(1:81, rep(1:27, each = 3))


### Construct a for loop for memory savings
tic()

for (i in seq_along(vec_DGM)){
  for (j in seq_along(vec_N_person)){
    for (k in seq_along(vec_WLE_rel)){
      
      # Print progress
      cat(paste0("DGM: ", vec_DGM[i]), 
          paste0("|| N_person:", vec_N_person[j]), 
          paste0("|| WLE_rel: ", vec_WLE_rel[k]),
          "\n")
      
      # Prepare parallel computation: Set the number of workers
      # parallelly::availableCores()
      # parallelly::availableWorkers()
      plan(multisession, workers = 8)
      
      # Load the saved dataframe
      file_name <- paste0("df", 
                          "_DGM-", vec_DGM[i], 
                          "_N-", vec_N_person[j], 
                          "_WLE_rel-", vec_WLE_rel[k]*100, ".rds")
      
      load_path <- file.path(data_dir2, 
                             "temp_folder", 
                             file_name)
      
      df_temp <- read_rds(load_path) %>%
        filter(beta_kind %in% c("beta_norm"))  # Extract only "beta_norm"
      
      # Prepare subset data for `pmap()`
      df_sub <- df_temp %>%
        dplyr::select(I, true_theta, vec_beta) %>%
        rename(
          N_item = I, 
          theta = true_theta, 
          beta = vec_beta) %>%
        mutate(N_person = map_dbl(.x = theta, .f = length)) %>%
        relocate(N_person, .before = N_item)
      
      # Simulate item response data (logit, prob, y)
      tic()
      
      list_sim_data <- future_pmap(
        .l = df_sub, 
        .f = gen_rasch_logit_prob_y,
        .options = furrr_options(seed = NULL, 
                                 chunk_size = 100, 
                                 scheduling = 5), 
        .progress = TRUE
      )
      
      toc()
      
      # Expand the item response lists
      df_y <- df_sub %>%
        mutate(
          logit = map(.x = list_sim_data, .f = 1), 
          prob = map(.x = list_sim_data, .f = 2), 
          y = map(.x = list_sim_data, .f = 3)
        )
      
      # Estimate test information
      df_y <- df_y %>%
        mutate(
          test_info_est = map_dbl(.x = logit, 
                                  .f = get_test_info_est)
      )
      
      # Attach to the original data
      df_temp2 <- df_temp %>% 
        dplyr::select(-N_person, -I)
      
      df_y_final <- bind_cols(df_temp2, df_y)
      
      
      # Save to the disk
      file_name <- paste0("df", 
                          "_DGM-", vec_DGM[i], 
                          "_N-", vec_N_person[j], 
                          "_WLE_rel-", vec_WLE_rel[k]*100, 
                          "_item response data.rds")
      
      save_path <- file.path(data_dir2, 
                             "temp_folder", 
                             file_name)
      
      write_rds(df_y_final, save_path)
    }
  }
}

toc()



###'#######################################################################'
###' 
###' Simulate item difficulty data 
###' 
###' via Monte Carlo approximation
###' 
###' (3) Fit the Rasch model and calculate WLE reliability
###'
###' - Loop over `3*5*5 = 75` files for memory savings
###' 
###' - Work only on the `beta_norm` for now to save time
###' 
###' - Sort out only the cases with `minimum WLE reliablities`
###' 
###'

### Prepare chopped looping for memory savings
vec_DGM <- c("Gaussian", "ALD", "Mixed")
vec_N_person <- c(20, 50, 100, 200, 500)
vec_WLE_rel <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# list_rep <- split(1:81, rep(1:27, each = 3))


### Prepare safe functions
safe_fit_rasch <- safely(fit_rasch) 
not_null <- negate(is_null)
safe_get_reliabilities <- safely(get_reliabilities)


### Construct a for loop for memory savings
tic()

for (i in seq_along(vec_DGM)){
  for (j in seq_along(vec_N_person)){
    for (k in seq_along(vec_WLE_rel)){
      
      # Print progress
      cat(paste0("DGM: ", vec_DGM[i]), 
          paste0("|| N_person:", vec_N_person[j]), 
          paste0("|| WLE_rel: ", vec_WLE_rel[k]),
          "\n")
      
      # Prepare parallel computation: Set the number of workers
      # parallelly::availableCores()
      # parallelly::availableWorkers()
      plan(multisession, workers = 8)
      
      # Load the saved dataframe
      file_name <- paste0("df", 
                          "_DGM-", vec_DGM[i], 
                          "_N-", vec_N_person[j], 
                          "_WLE_rel-", vec_WLE_rel[k]*100, 
                          "_item response data.rds")
      
      load_path <- file.path(data_dir2, 
                             "temp_folder", 
                             file_name)
      
      df_temp <- read_rds(load_path) %>%
        filter(N_item >= 5)  ## Remove cases with N_item < 5
      
      # Fit the Rasch model
      tic()
      
      plan(multisession, workers = 12)
      
      df_fit <- df_temp %>%
        # slice_head(n = 2000) %>%
        mutate(
          mod_fit = future_map(.x = y, 
                               .f = safe_fit_rasch, 
                               .options = furrr_options(seed = NULL, 
                                                        chunk_size = 200, 
                                                        scheduling = 5), 
                               .progress = TRUE)
        )
      
      toc()
      
      # Split result and error
      df_fit <- df_fit %>%
        mutate(
          error = map(.x = mod_fit, .f = "error"), 
          result = map(.x = mod_fit, .f = "result"), 
          error_lgl = map_lgl(.x = error, .f = not_null)
        ) %>%
        dplyr::select(-mod_fit)
      
      # Calculate reliabilities 
      tic()
      
      plan(multisession, workers = 8)
      
      df_rel <- df_fit %>%
        # slice_head(n = 2000) %>%
        mutate(
          rel_est_temp = future_map(.x = result, 
                                    .f = safe_get_reliabilities, 
                                    .options = furrr_options(seed = NULL, 
                                                             chunk_size = 500, 
                                                             scheduling = 1), 
                                    .progress = TRUE), 
          rel_est = map(.x = rel_est_temp, .f = "result"), 
          WLE_rel_est = map(.x = rel_est, .f = "WLE_rel") %>%
            map_dbl(.f = ~if_else(is_null(.), NA_real_, .)), 
          EAP_rel_est = map(.x = rel_est, .f = "EAP_rel") %>%
            map_dbl(.f = ~if_else(is_null(.), NA_real_, .))
        ) %>%
        dplyr::select(-rel_est_temp, -rel_est)
      
      toc()
      
      # Calculate biases
      df_rel <- df_rel %>%
        mutate(
          WLE_rel_bias = WLE_rel_est - WLE_rel, 
          EAP_rel_bias = EAP_rel_est - EAP_rel,
          test_info_bias = test_info_est - test_info
        )
      
      #' Filter out results based on WLE reliability biases
      #' reduce rows
      df_rel_min <- df_rel %>%
        # drop errors
        filter(error_lgl == FALSE) %>% 
        group_by(sim_cond, rep) %>%
        mutate(
          min_WLE_bias = min(abs(WLE_rel_bias))
        ) %>%
        ungroup() %>%
        filter(min_WLE_bias == abs(WLE_rel_bias))
      
      
      # Save to the disk
      file_name <- paste0("df", 
                          "_DGM-", vec_DGM[i], 
                          "_N-", vec_N_person[j], 
                          "_WLE_rel-", vec_WLE_rel[k]*100, 
                          "_final data with min_WLE_rel.rds")
      
      save_path <- file.path(data_dir2, 
                             "temp_folder", 
                             file_name)
      
      write_rds(df_rel_min, save_path)
    }
  }
}

toc()


