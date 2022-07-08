
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Results
###' 
###' Task: Shiny App Development
###'       `Best choice` 
###'       
###' Data: Performance evaluator estimates (MSEL, MSELP, ISEL, etc)
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
library(rlang)
library(broom)
library(ggeffects)
library(glue)

library(sandwich)
library(lmtest)
library(miceadds)
library(multcomp)
library(clubSandwich)
library(margins)
library(effects)
library(cowplot)
library(gtools)


### Call custom functions
list.files(file.path(work_dir, "functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################
###'
###' Load the performance evaluator files
###'
###'

load_path_wide <- file.path(data_dir2, "df_loss_estimates_WIDE_collected.rds")
load_path_long <- file.path(data_dir2, "df_loss_estimates_LONG_collected.rds")
load_path_sum <- file.path(data_dir2, "df_loss_estimates_SUMMARY_collected.rds")

df_wide <- read_rds(load_path_wide)
df_long <- read_rds(load_path_long)
df_sum <- read_rds(load_path_sum)



###'#######################################################################
###'
###' `best_choice_map()`
###' 
###' - Best choice heatmap according to the five data factors
###' - Indicate the location of the selected data factors on the heatmap
###'
###'

best_choice_map <- function(df_sum = df_sum,
                            loss_est_fix = "MSEL", 
                            param_fix = "theta", 
                            DGM_fix = "Mixed",
                            N_person_fix = "N = 100",
                            WLE_rel_fix = "WLE rel. = 0.6"){
  
  # Filter only relevant loss estimates
  df_sum2 <- df_sum %>%
    filter(loss_est == loss_est_fix)
  
  # Sort out only the best choice raws
  vec_groups <- c("param", "DGM", "N_person", "WLE_rel")
  
  df_best <- df_sum2 %>%
    group_by(
      across(.cols = all_of(vec_groups))
    ) %>%
    arrange(
      across(.cols = all_of(c(vec_groups, "mean_raw")))
    ) %>%
    mutate(min_loss = min(mean_raw)) %>%
    ungroup() %>%
    filter(mean_raw == min_loss) 
    
  # Change factor labels into more brief version
  levels(df_best$WLE_rel) <- levels(df_best$WLE_rel) %>%
    str_remove("WLE rel. = ")

  levels(df_best$N_person) <- levels(df_best$N_person) %>%
    str_remove("N = ")
  
  # Labels
  # labs(y = TeX("$\\sigma_\\tau$"),
  #      x = TeX("$\\bar{n}_j$"),
  #      color = "Model",
  #      pch = "Summary\nmethod")
  
  subtitle_text <- glue("Selected: ", 
                        "parameter = {param_fix}, ", 
                        "True G = {DGM_fix}, ", 
                        "{N_person_fix}, ", 
                        "{WLE_rel_fix}")
  
  labels  <- labs(
    title = glue("{loss_est_fix}: Best-choice Map"),
    subtitle = subtitle_text, 
    x = "WLE reliability", 
    y = "Number of Persons (N)", 
    caption = "Column: Latent parameter; Row: True G dist."
  )
  
  # Selected label text
  ann_text <- df_best %>%
    filter(
      .data[["param"]] == param_fix, 
      .data[["DGM"]] == DGM_fix, 
      .data[["N_person"]] == N_person_fix %>% str_remove("N = "), 
      .data[["WLE_rel"]] == WLE_rel_fix %>% str_remove("WLE rel. = ")
    ) %>%
    dplyr::select(all_of(vec_groups))
  
  
  # Plot heatmap
  df_best %>%
    ggplot() +
    geom_point(aes(x = WLE_rel, y= N_person, color = model, pch = sum_method),
               size = 6) +
    scale_shape_discrete(drop = F) +
    facet_grid(
      rows = vars(DGM), 
      cols = vars(param)
    ) +
    theme_bw() + labels + 
    theme(legend.box = "horizontal") + 
    geom_text(data = ann_text,
              mapping = aes(x = WLE_rel,
                            y = N_person,
                            label = "Selected"),
              vjust = -0.6)
}

### Test the function
best_choice_map(df_sum = df_sum, 
                param_fix = "theta", 
                DGM_fix = "Mixed",
                N_person_fix = "N = 100",
                WLE_rel_fix = "WLE rel. = 0.9", 
                loss_est_fix = "IAEL")



###'#######################################################################
###'
###' `best_choice_table()`
###'
###' - Best choice descriptive statistics (mean comparisons)
###' 
###'

best_choice_table <- function(df_sum = df_sum, 
                              # loss_est_fix = "MSEL", 
                              param_fix = "theta", 
                              DGM_fix = "Mixed",
                              N_person_fix = "N = 100",
                              WLE_rel_fix = "WLE rel. = 0.6") {
  
  # Filter/Subset the data
  df_sub <- df_sum %>%
    filter(
      # loss_est == loss_est_fix, 
      param == param_fix, 
      DGM == DGM_fix, 
      N_person == N_person_fix, 
      WLE_rel == WLE_rel_fix
    ) %>%
    arrange(loss_est, mean_raw)
  
  # Calculate proportional benefit
  df_sub2 <- df_sub %>%
    group_by(loss_est) %>%
    mutate(
      # take means
      mean_raw_min = min(mean_raw), 
      mean_log_min = min(mean_log), 
      
      # proportional change
      delta_raw = mean_raw - mean_raw_min, 
      pct_delta_raw = 100*(delta_raw/mean_raw_min), 
      
      delta_log = mean_log - mean_log_min, 
      pct_delta_log_approx = 100*delta_log
    )
  
  # Return output
  df_sub2
}


### Test
df_best <- best_choice_table(
  df_sum = df_sum,
  # loss_est_fix = "ISEL", 
  param_fix = "theta", 
  DGM_fix = "Mixed",
  N_person_fix = "N = 100",
  WLE_rel_fix = "WLE rel. = 0.6")

df_temp <- df_best %>%
  mutate(across(where(is.numeric), ~round(., 3)))

View(df_temp)



###'#######################################################################
###'
###' `best_choice_cond_pred()`
###'
###' - Best choice conditional predictions model 
###'   (with cluster robust standard errors)
###'
###'

best_choice_cond_pred <- function(df_long = df_long, 
                                  loss_est_fix = "MSEL", 
                                  param_fix = "theta", 
                                  DGM_fix = "Mixed",
                                  N_person_fix = "N = 100",
                                  WLE_rel_fix = "WLE rel. = 0.6"){
  
  # Filter/Subset the data
  data <- df_long %>%
    filter(
      loss_est == loss_est_fix, 
      param == param_fix, 
      DGM == DGM_fix, 
      N_person == N_person_fix, 
      WLE_rel == WLE_rel_fix
    ) %>%
    arrange(cluster_ID, loss_est, model, sum_method)
  
  # Extract cluster ID
  vec_cluster_ID <- data$cluster_ID
  
  # Fit the model 
  lm_fit <- lm(formula = "value ~ (model + sum_method)^2", 
               data = data)
  
  # Generate conditional predictions
  df_cond_pred <- ggpredict(model = lm_fit, 
                            terms = c("model", "sum_method"),
                            vcov.fun = "vcovCR", 
                            vcov.type = "CR0", 
                            vcov.args = list(cluster = vec_cluster_ID)) %>%
    tibble() %>%
    mutate(loss_est = loss_est_fix) %>%
    dplyr::select(loss_est, everything())
  
  
  # Extract label information
  subtitle_lab <- paste(
    paste0("parameter: ", param_fix), 
    paste0("True G: ", DGM_fix), 
    N_person_fix, 
    WLE_rel_fix,
    collapse = "", sep = ", "
  )
  
  labels <- 
    labs(x = "Model for G", 
         y = paste0("Predicted (", loss_est_fix, ")"), 
         title = paste0("Predicted ",  loss_est_fix," by meta−model regression"),
         subtitle = subtitle_lab, 
         caption = NULL)
  
  # Create a plotdd
  # Generate the main plot (without panel)
  p <- ggplot(data = df_cond_pred, 
              aes(x = x, y = predicted, 
                  group = group, color = group, shape = group)) +
    
    geom_point(aes(y = predicted), size = 3, position = position_dodge(width = 0.4)) + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0.4), width = 0.2) + 
    geom_line(aes(y = predicted), position = position_dodge(width = 0.4)) + 
    # geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
    
    geom_text(aes(label = sprintf("%#.3f", round(predicted, 3)) %>% 
                    str_replace("0.", ".")),
              hjust = 1.5, color = "black", size = 3, 
              position = position_dodge(width = 0.4)) +
    
    scale_x_discrete(expand = expansion(mult = 0.4)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) + 
    
    scale_color_manual(values = color_palette[seq(unique(df_cond_pred$group))]) + 
    scale_shape_manual(values = shape_palette[seq(unique(df_cond_pred$group))]) + 
    
    theme_trend + 
    theme(strip.background = element_rect(fill = "gray100")) + 
    labels
  
  # Return list
  list(loss_est_fix, subtitle_lab, data, df_cond_pred, p)
}

### Test
best_choice_cond_pred(df_long = df_long, 
                      loss_est_fix = "KS_dist", 
                      param_fix = "theta", 
                      DGM_fix = "Mixed",
                      N_person_fix = "N = 200",
                      WLE_rel_fix = "WLE rel. = 0.7")



###'#######################################################################
###'
###' `best_choice_wrapper()`
###'
###' Define a wrapper function to get the best choice results
###'
###' - vec_param has to be a named vector with an order
###'
###'

best_choice_wrapper <- function(vec_param, df_sum, df_long){
  
  # Generate the Best choice summary
  df_best <- best_choice_table(
    df_sum = df_sum,
    # loss_est_fix = vec_param["loss_est_fix"],
    param_fix = vec_param["param_fix"], 
    DGM_fix = vec_param["DGM_fix"],
    N_person_fix = vec_param["N_person_fix"],
    WLE_rel_fix = vec_param["WLE_rel_fix"]
  )
  
  # Tidy up as a table
  df_best_table <- df_best %>%
    dplyr::select(loss_est, model, sum_method, 
                  mean_raw, delta_raw, pct_delta_raw) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  # Collect conditional predictions across performance evaluators
  vec_loss_est <- df_best_table$loss_est %>% 
    unique()
  
  list_collect <- list()
  
  for (i in seq_along(vec_loss_est)){
    
    # Get the list of conditional predictions 
    list_best_cond_pred <- 
      best_choice_cond_pred(df_long = df_long, 
                            loss_est_fix = vec_loss_est[i],
                            param_fix = vec_param["param_fix"], 
                            DGM_fix = vec_param["DGM_fix"],
                            N_person_fix = vec_param["N_person_fix"],
                            WLE_rel_fix = vec_param["WLE_rel_fix"])
    
    # Append the df_best and df_best_table
    list_best_cond_pred[[6]] <- df_best %>%
      filter(loss_est == vec_loss_est[i])
    
    list_best_cond_pred[[7]] <- df_best_table %>%
      filter(loss_est == vec_loss_est[i])
    
    list_collect[[i]] <- list_best_cond_pred
  }
  
  # Return output
  list_collect
}


### Test
vec_param <- c(
  loss_est_fix = "KS_dist", 
  param_fix = "theta", 
  DGM_fix = "Mixed",
  N_person_fix = "N = 200",
  WLE_rel_fix = "WLE rel. = 0.7"
)

list_temp <- best_choice_wrapper(vec_param, df_sum, df_long)

list_temp

map(list_temp, 7)

