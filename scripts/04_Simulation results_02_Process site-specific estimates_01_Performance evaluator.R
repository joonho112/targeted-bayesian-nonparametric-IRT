
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Results
###' 
###' Task: Process site-specific estiamtes
###'        `Loss estimates - performance evaluators`
###'       
###' Data: Site-specific estimates
###' 
###' Date: 
###' - 2022-07-05
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
###' Load / Collect all site-specific estimates
###'
###'

load_path <- file.path(data_dir2, 
                       "df_site_est", 
                       "df_site_est_collected.rds")

df_site_est <- read_rds(load_path)

object.size(df_site_est) %>% format("MB")



###'#######################################################################'
###'
###' Apply performance evaluator function and
###' get performance summaries
###' 
###' `get_losses()`
###' 
###'  - `MSEL`
###'  - `MSELR`
###'  - `MSELP`
###'  - `ISEL`
###'  - `ISEL_hh`
###'  - `KS_dist`
###'
###'

get_losses <- function(df_est) {
  
  # (1) Detect the parameter name
  param <- names(df_est) %>%
    str_detect("_true") %>%
    names(df_est)[.] %>%
    str_remove("_true")
  
  # (2) Generate a nested data
  df_nest <- df_est %>%
    dplyr::select(contains(c("true", "pm", "cb", "gr"))) %>%
    pivot_longer(
      cols = contains(c("pm", "cb", "gr")), 
      names_to = "sum_method", names_prefix = paste0(param, "_"), 
      values_to = "estimate"
    ) %>%
    mutate(
      sum_method = factor(sum_method, 
                          levels = c("pm", "cb", "gr"), 
                          labels = c("PM", "CB", "GR"))
    ) %>%
    set_names(c("true", "sum_method", "estimate")) %>%
    arrange(sum_method) %>%
    group_by(sum_method) %>%
    nest() %>%
    ungroup()
  
  
  # (3) Apply loss functions to each nested dataset
  df_loss <- df_nest %>%
    mutate(
      MSEL = map_dbl(.x = data, .f = ~MSEL(.x$estimate, .x$true)), 
      MSELR = map_dbl(.x = data, .f = ~MSELR(.x$estimate, .x$true)),
      MSELP = map_dbl(.x = data, .f = ~MSELP(.x$estimate, .x$true)),
      ISEL = map_dbl(.x = data, .f = ~ISEL(.x$estimate, .x$true)),
      IAEL = map_dbl(.x = data, .f = ~IAEL(.x$estimate, .x$true)),
      KS_dist = map_dbl(data, ~KS_dist(.x$estimate, .x$true))
    ) %>%
    select(-data)
  
  return(df_loss)
}


###' Test the function
###' Calculate performance indicators: MSEL, MSELP, ISEL etc 
df_loss_theta <- df_site_est %>%
  slice_head(n = 10) %>%  # for testing
  mutate(
    loss_theta = map(.x = df_theta, 
                     .f = get_losses)
  ) %>%
  select(-(theta:y), -(file_name:file_path), -(df_theta:df_hyper)) %>%
  unnest(loss_theta)

df_loss_beta <- df_site_est %>%
  slice_head(n = 10) %>%  # for testing
  mutate(
    loss_beta = map(.x = df_beta, 
                    .f = get_losses)
  ) %>%
  select(-(theta:y), -(file_name:file_path), -(df_theta:df_hyper)) %>%
  unnest(loss_beta)



###'#######################################################################'
###'
###' Apply the defined function for performance evaluation
###'
###'

### Prepare parallel computation: Set the number of workers
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 19)


### Calculate performance indicators: (1) theta
tic()

df_loss_theta <- df_site_est %>%
  # slice_head(n = 10) %>%  # for testing
  mutate(
    loss_theta = future_map(.x = df_theta, 
                            .f = get_losses,
                            .options = furrr_options(seed = NULL,
                                                     chunk_size = 100,
                                                     scheduling = 2),
                            .progress = TRUE)
  ) %>%
  select(-(theta:y), -(file_name:file_path), -(df_theta:df_hyper)) %>%
  unnest(loss_theta) %>%
  mutate(param = "theta")

toc()


### Calculate performance indicators: (2) beta
tic()

df_loss_beta <- df_site_est %>%
  # slice_head(n = 10) %>%  # for testing
  mutate(
    loss_beta = future_map(.x = df_beta, 
                           .f = get_losses,
                           .options = furrr_options(seed = NULL,
                                                    chunk_size = 100,
                                                    scheduling = 2),
                           .progress = TRUE)
  ) %>%
  select(-(theta:y), -(file_name:file_path), -(df_theta:df_hyper)) %>%
  unnest(loss_beta) %>%
  mutate(param = "beta")

toc()


### Combine theta and beta results
df_loss <- bind_rows(df_loss_theta, df_loss_beta) %>%
  mutate(
    param = factor(param, levels = c("theta", "beta"))
  ) 


### Arrange the resulting dataset
df_loss2 <- df_loss %>% 
  select(param, sim_cond, cond_name, DGM, N_person, true_var, WLE_rel, N_item, rep, 
         model, sum_method, everything()) %>%
  arrange(param, cond_name, model, sum_method) %>%
  relocate(MSEL:KS_dist, .after = "sum_method")

object.size(df_loss2) %>% format("MB")



###'#######################################################################'
###'
###' Generate a cluster ID
###'
###' - We want to apply `cluster robust standard errors`
###'   because multiple models were fitted to the SAME data
###'
###' - Cluster ID = Pre-generated data ID
###' 
###' - Group by: 
###'   (DGM, N_person, WLE_rel) + rep
###'   sim_cond + rep
###'   
###' - 9 rows (3 model by 3 posterior summary) per cluster
###' - 18 rows for theta + beta
###'
###'

df_loss3 <- df_loss2 %>%
  group_by(sim_cond, rep) %>%
  mutate(cluster_ID = cur_group_id()) %>%
  relocate(cluster_ID, .after = sim_cond) %>%
  arrange(param, sim_cond, rep, model, sum_method) %>%
  ungroup() %>%
  select(-beta_kind)



###'#######################################################################
###'
###' Save the loss estimates - `WIDE` format data
###'
###'

save_path <- file.path(data_dir2, 
                       "df_site_est", 
                       "df_loss_estimates_WIDE_collected.rds")

write_rds(df_loss3, save_path)



###'#######################################################################
###'
###' Generate a reference table
###' 
###' 

vec_var <- c("DGM", "N_person", "WLE_rel",  
             "model", "sum_method")

list_levels <- list(
  lev_DGM = c("Gaussian", "ALD", "Mixed"),
  lev_N_person = c(20, 50, 100, 200, 500), 
  lev_WLE_rel = c(0.50, 0.60, 0.70, 0.80, 0.90), 
  lev_model = c("Gaussian", "DP-inform", "DP-diffuse"), 
  lev_sum_method = c("PM", "CB", "GR")
)

list_labels <- list(
  lab_DGM = list_levels[["lev_DGM"]],
  lab_N_person = paste0("N = ", list_levels[["lev_N_person"]]), 
  lab_WLE_rel = paste0("WLE rel. = ", list_levels[["lev_WLE_rel"]]), 
  lab_model = list_levels[["lev_model"]], 
  lab_sum_method = list_levels[["lev_sum_method"]] 
)

tbl_ref <- tibble(vec_var, list_levels, list_labels) %>%
  set_names(c("factor", "level", "label"))



###'#######################################################################
###'
###' Convert to a `LONG` format data
###'
###'

### Load the wide file
load_path <- file.path(data_dir2, 
                       "df_site_est", 
                       "df_loss_estimates_WIDE_collected.rds")

df_loss_wide <- read_rds(load_path)

object.size(df_loss_wide) %>% format("MB")


### Convert to long data format 
df_loss_wide %>% 
  ungroup() %>%
  dplyr::select(MSEL:KS_dist) %>%
  names() -> vec_loss_est

df_loss_long <- df_loss_wide %>%
  ungroup() %>%
  select(-(test_info_est:min_WLE_bias)) %>%
  pivot_longer(cols = MSEL:KS_dist, 
               names_to = "loss_est", 
               values_to = "value") %>%
  mutate(loss_est = factor(loss_est, 
                           levels = vec_loss_est))

View(df_loss_long)


### Convert to factors
df_loss_long2 <- df_loss_long %>%
  mutate(

    N_person = factor(
      N_person, 
      levels = list_levels[["lev_N_person"]], 
      labels = list_labels[["lab_N_person"]]
    ), 
    
    WLE_rel = factor(
      WLE_rel, 
      levels = list_levels[["lev_WLE_rel"]], 
      labels = list_labels[["lab_WLE_rel"]]
    ), 
    
    N_item = paste0("I = ", N_item)
    
  ) %>%
  arrange(param, cluster_ID, loss_est, model, sum_method)



###'#######################################################################
###'
###' Save the loss estimates - `LONG` format data
###'
###'

save_path <- file.path(data_dir2, 
                       "df_site_est", 
                       "df_loss_estimates_LONG_collected.rds")

object.size(df_loss_long2) %>% format("MB")

write_rds(df_loss_long2, save_path)

