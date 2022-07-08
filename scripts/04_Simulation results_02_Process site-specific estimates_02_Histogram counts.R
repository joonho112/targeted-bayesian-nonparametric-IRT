
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Results
###' 
###' Task: Process site-specific estiamtes
###'        `Histogram Counts`
###'       
###' Data: Site-specific estimates
###' 
###' Date: 
###' - 2022-07-06
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
###' Load / Collect all site-specific estimates
###'
###'

# ### Load only the collected data
# load_path <- file.path(data_dir2, 
#                        "df_site_est_collected.rds")
# 
# df_site_est <- read_rds(load_path)
# 
# object.size(df_site_est) %>% format("MB")


### Collect all dataframes
save_path <- file.path(data_dir2, "df_site_est")

tic()

df_temp <- list.files(save_path) %>%
  tibble() %>%
  set_names(c("file_name")) %>%
  mutate(file_path = file.path(save_path, file_name)) %>%
  mutate(
    data = future_map(.x = file_path, 
                      .f = read_rds, 
                      .options = furrr_options(seed = NULL,
                                               chunk_size = 1,
                                               scheduling = 1),
                      .progress = TRUE)
  )

toc()


### Bind rows
df_site_est <- bind_rows(df_temp$data) %>%
  arrange(model, sim_cond, rep)

object.size(df_site_est) %>% format("GB")



###'#######################################################################'
###'
###' Get histogram counts
###' 
###' `get_hist_count()`
###'
###' - Approximate super-population distribution
###' 
###' 

get_hist_count <- function(df_est){
  
  # (1) Define a function to count just for one site-specific parameter vector
  one_hist_count <- function(vec_est){
    
    # Prepare bins - by default, 50 bins for (-6*SD, +6*SD)
    cutpoints <- seq(from = -6, to = 6, length.out = 51)
    
    # Generate a tibble containing only one count vector
    df_count <- cut(vec_est, breaks = cutpoints) %>%
      table() %>% 
      as.data.frame() %>%
      rownames_to_column() %>%
      dplyr::select(-rowname) %>%
      set_names(c("bin", "count")) %>%
      tibble()
    
    return(df_count)
  }
  
  # (2) Apply the function to tau_j vectors
  est_lev <- c("true", "pm", "cb", "gr")
  
  df_vec_est <- df_est %>% 
    dplyr::select(contains(est_lev)) 
  
  df_hist <- map_dfr(.x = df_vec_est,  
                     .f = one_hist_count, 
                     .id = "estimator") 
  
  return(df_hist)
}



###'#######################################################################
###'
###' Calculate histogram counts
###' 
###' Apply the function - `get_hist_count()`
###' 
###' 

### Calculate histogram counts: (1) theta
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 13)

tic()

df_his_theta <- df_site_est %>%
  # slice_head(n = 10) %>%
  select(sim_cond, cond_name, DGM, N_person, true_var, WLE_rel, N_item, rep, model, 
         df_theta) %>%
  
  mutate(
    
    hist = future_map(.x = df_theta,
                      .f = get_hist_count, 
                      .options = furrr_options(seed = NULL,
                                               chunk_size = 1000,
                                               scheduling = 1),
                      .progress = TRUE)
    
  ) %>%
  dplyr::select(-df_theta)

toc()


### Calculate histogram counts: (2) beta
parallelly::availableCores()
parallelly::availableWorkers()
plan(multisession, workers = 13)

tic()

df_his_beta <- df_site_est %>%
  # slice_head(n = 10) %>%
  select(sim_cond, cond_name, DGM, N_person, true_var, WLE_rel, N_item, rep, model, 
         df_beta) %>%
  
  mutate(
    
    hist = future_map(.x = df_beta,
                      .f = get_hist_count, 
                      .options = furrr_options(seed = NULL,
                                               chunk_size = 1000,
                                               scheduling = 1),
                      .progress = TRUE)
    
  ) %>%
  dplyr::select(-df_beta)

toc()



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
###' Tidy up the calculated histogram counts - (1) `theta`
###'
###'

###' (1) Unnest and split estimator into `param` and `sum_method`
tic()

df_his_unnest <- df_his_theta %>%
  unnest(hist) %>%
  # slice_head(n = 100000) %>%
  separate(col = estimator, into = c("param", "sum_method"), sep = "_")

toc()


###' (2) Calculate the group sums
df_his_sum <- df_his_unnest %>%
  group_by(param, DGM, N_person, WLE_rel, model, sum_method, bin) %>%
  summarize(
    N_sum = sum(count),
    N_rep = n_distinct(rep)
  ) %>% 
  ungroup()


###' (3) Calculate the scaled EDF
df_his_scaled <- df_his_sum %>%
  group_by(param, DGM, N_person, WLE_rel, model, sum_method) %>%
  mutate(
    N_total = N_person*N_rep, 
    density = N_sum/N_total
  )


###' (4) Split Bin into start and end points
df_hist <- df_his_scaled$bin %>%
  str_split_fixed(",", n = 2) %>%
  data.frame() %>%
  mutate(start = str_sub(X1, start = 2) %>% as.numeric(), 
         end = str_sub(X2, end = -2) %>% as.numeric(), 
         middle = (start + end)/2) %>%
  dplyr::select(-X1, -X2) %>%
  cbind.data.frame(df_his_scaled) %>%
  relocate(start, end, middle, .after = bin) %>%
  tibble() 


###' (5) Generate a density column for `true` values
group_key <- c("DGM", "N_person", "WLE_rel", "model", 
               "bin")

df_true_dens <- df_hist %>% 
  filter(sum_method == "true") %>%
  dplyr::select(all_of(group_key), density) %>%
  rename(true_dens = density)

df_join <- df_hist %>%
  full_join_track(df_true_dens, by = group_key) 


### (6) Assign factor levels & Arrange
df_hist_theta <- df_join %>%
  mutate(
    
    DGM = factor(
      DGM, 
      levels = list_levels[["lev_DGM"]]
    ), 
    
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
    
    model = factor(
      model, 
      levels = list_levels[["lev_model"]]
    ), 
    
    sum_method = factor(
      sum_method, 
      levels = c("true", "pm", "cb", "gr"),  
      labels = c("true", list_levels[["lev_sum_method"]])
    )
  ) %>%
  arrange(DGM, N_person, WLE_rel, model, sum_method, bin) %>%
  select(param:middle, N_sum, density, true_dens, N_rep, N_total)




###'#######################################################################
###'
###' Tidy up the calculated histogram counts - (2) `beta`
###'
###'

###' (1) Unnest and split estimator into `param` and `sum_method`
tic()

df_his_unnest <- df_his_beta %>%
  unnest(hist) %>%
  # slice_head(n = 100000) %>%
  separate(col = estimator, into = c("param", "sum_method"), sep = "_")

toc()


###' (2) Calculate the group sums
df_his_sum <- df_his_unnest %>%
  group_by(param, DGM, N_person, WLE_rel, model, sum_method, bin) %>%
  summarize(
    N_sum = sum(count),
    N_rep = n_distinct(rep), 
    N_total = sum(N_item)
  ) %>% 
  ungroup()


###' (3) Calculate the scaled EDF
df_his_scaled <- df_his_sum %>%
  group_by(param, DGM, N_person, WLE_rel, model, sum_method) %>%
  mutate(
    density = N_sum/N_total
  )


###' (4) Split Bin into start and end points
df_hist <- df_his_scaled$bin %>%
  str_split_fixed(",", n = 2) %>%
  data.frame() %>%
  mutate(start = str_sub(X1, start = 2) %>% as.numeric(), 
         end = str_sub(X2, end = -2) %>% as.numeric(), 
         middle = (start + end)/2) %>%
  dplyr::select(-X1, -X2) %>%
  cbind.data.frame(df_his_scaled) %>%
  relocate(start, end, middle, .after = bin) %>%
  tibble() 


###' (5) Generate a density column for `true` values
group_key <- c("DGM", "N_person", "WLE_rel", "model", 
               "bin")

df_true_dens <- df_hist %>% 
  filter(sum_method == "true") %>%
  dplyr::select(all_of(group_key), density) %>%
  rename(true_dens = density)

df_join <- df_hist %>%
  full_join_track(df_true_dens, by = group_key) 


### (6) Assign factor levels & Arrange
df_hist_beta <- df_join %>%
  mutate(
    
    DGM = factor(
      DGM, 
      levels = list_levels[["lev_DGM"]]
    ), 
    
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
    
    model = factor(
      model, 
      levels = list_levels[["lev_model"]]
    ), 
    
    sum_method = factor(
      sum_method, 
      levels = c("true", "pm", "cb", "gr"),  
      labels = c("true", list_levels[["lev_sum_method"]])
    )
  ) %>%
  arrange(DGM, N_person, WLE_rel, model, sum_method, bin) %>%
  select(param:middle, N_sum, density, true_dens, N_rep, N_total)



###'#######################################################################
###'
###' Bind rows and save as .rds file
###'
###'

df_hist_theta
df_hist_beta

df_hist <- bind_rows(df_hist_theta, df_hist_beta)

save_path <- file.path(data_dir2, "df_histogram_collected.rds")

write_rds(df_hist, save_path)



###'#######################################################################
###'
###' Add "Super-population" true density
###' 
###' based on theoretical Data-Generating distribution
###'
###' - `theta`: based on DGM
###' - `beta`: regardless of DGM
###'
###'

### Load the previous version of histogram counts
load_path <- file.path(data_dir2, "df_histogram_collected.rds")
df_hist <- read_rds(load_path)

View(df_hist)

df_hist %>%
  count(param, DGM)


### Simulate theoretical densities with extremely large sized sampling
N_large <- 100000

vec_beta <- rnorm(n = N_large, mean = 0, sd = 1)

vec_theta_Gaussian <- gen_vec_theta(DGM = "Gaussian", N = N_large)

vec_theta_ALD <- gen_vec_theta(DGM = "ALD", N = N_large)

vec_theta_Mixed <- gen_vec_theta(DGM = "Mixed", N = N_large)

plot(density(vec_beta))
plot(density(vec_theta_Gaussian))
plot(density(vec_theta_ALD))
plot(density(vec_theta_Mixed))


###' Generate histogram counts based on the theoretical true densities
###' Define a function to count just for one site-specific parameter vector
one_hist_true_dens <- function(vec_est, N_large = 100000){
  
  # Prepare bins - by default, 50 bins for (-6*SD, +6*SD)
  cutpoints <- seq(from = -6, to = 6, length.out = 51)
  
  # Generate a tibble containing only one count vector
  df_true_dens <- cut(vec_est, breaks = cutpoints) %>%
    table() %>% 
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::select(-rowname) %>%
    set_names(c("bin", "count")) %>%
    tibble() %>%
    mutate(true_dens = count/N_large)
  
  return(df_true_dens)
}


### Prepare a theoretical true density tibbles to merge into original dataset
df_hist

# (1) beta - overall
df_true_beta <- one_hist_true_dens(vec_beta, N_large) %>%
  mutate(param = "beta") %>%
  dplyr::select(-count) %>%
  dplyr::select(param, bin, true_dens)


# (2) theta - Gaussian
df_true_theta_Gaussian <- one_hist_true_dens(vec_theta_Gaussian, N_large) %>%
  mutate(param = "theta", DGM = "Gaussian")

# (3) theta - ALD
df_true_theta_ALD <- one_hist_true_dens(vec_theta_ALD, N_large) %>%
  mutate(param = "theta", DGM = "ALD")

# (4) theta - Mixed
df_true_theta_Mixed <- one_hist_true_dens(vec_theta_Mixed, N_large) %>%
  mutate(param = "theta", DGM = "Mixed")

# (5) Combine into one theta tibble
df_true_theta <- bind_rows(df_true_theta_Gaussian, 
                           df_true_theta_ALD, 
                           df_true_theta_Mixed) %>%
  dplyr::select(-count) %>%
  relocate(param, DGM, .before = "bin") %>%
  mutate(DGM = factor(DGM, levels = c("Gaussian", "ALD", "Mixed")))


### Merge into the original df_hist dataframe
df_hist

df_hist_theta <- df_hist %>%
  dplyr::select(-true_dens) %>%
  filter(param == "theta") %>%
  full_join_track(df_true_theta, by = c("param", "DGM", "bin"), .merge = FALSE)

df_hist_beta <- df_hist %>%
  dplyr::select(-true_dens) %>%
  filter(param == "beta") %>% 
  full_join_track(df_true_beta, by = c("param", "bin"), .merge = FALSE)

df_hist2 <- bind_rows(df_hist_theta, df_hist_beta) %>%
  mutate(param = factor(param, levels = c("theta", "beta"))) %>%
  relocate(true_dens, .after = "density")

View(df_hist2)


### Save the new dataframe
save_path <- file.path(data_dir2, "df_histogram_collected_new_true_dens.rds")
write_rds(df_hist2, save_path)
