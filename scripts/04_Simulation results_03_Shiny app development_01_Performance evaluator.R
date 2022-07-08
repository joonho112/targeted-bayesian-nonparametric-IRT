
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Results
###' 
###' Task: Shiny App Development
###'       `Performance evaluators` 
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
###' Load the performance evaluator files
###'
###'

load_path_wide <- file.path(data_dir2, "df_loss_estimates_WIDE_collected.rds")
load_path_long <- file.path(data_dir2, "df_loss_estimates_LONG_collected.rds")

df_wide <- read_rds(load_path_wide)
df_long <- read_rds(load_path_long) %>%
  mutate(DGM = factor(DGM, levels = c("Gaussian", "ALD", "Mixed")))



###'#######################################################################
###'
###' Define the total SEVEN factors and their levels
###'
###'

vec_7_factors <- c("param", "DGM", "N_person", "WLE_rel", 
                 "model", "sum_method", "loss_est")

list_levels <- df_long %>% 
  dplyr::select(all_of(vec_7_factors)) %>%
  map(.f = levels)



###'#######################################################################
###'
###' `fix_3_factor()`
###' 
###' A function to fix three factors
###' 
###' -> .data pronoun version
###'
###'

fix_3_factors <- function(
    df = df_long, 
    var1 = "loss_est", lev1 = "MSEL",
    var2 = "param", lev2 = "theta",
    var3 = "N_person", lev3 = "N = 100"){
  
  # Create a vector of fixed factors
  vec_fixed <- c(var1, var2, var3)
  
  # Create a list of fixed levels
  list_fixed <- list(
    paste0(var1, ": ", lev1), 
    paste0(var2, ": ", lev2), 
    paste0(var3, ": ", lev3)
  ) %>%
    unlist() %>%
    paste(collapse = ", ")
  
  # Filter the original dataframe
  df_sub <- df %>%
    filter(.data[[var1]] == lev1) %>%
    filter(.data[[var2]] == lev2) %>%
    filter(.data[[var3]] == lev3) 
  
  # Summarize log means
  vec_7_factors <- c("param", "DGM", "N_person", "WLE_rel", 
                     "model", "sum_method", "loss_est")
  
  df_sub2 <- df_sub %>%
    group_by(
      across(.cols = all_of(vec_7_factors))
    ) %>%
    mutate(
      log_value = if_else(value == 0, 0, log(value*1000))
    ) %>%
    summarize(
      mean_log = mean(log_value), 
      .groups = 'drop'
    ) %>%
    dplyr::select(-all_of(vec_fixed))
  
  # Return the resulting values
  list(vec_fixed, list_fixed, df_sub2)
}

### Test the function
list_sub <- fix_3_factors(
  df_long, 
  var1 = "loss_est", lev1 = "MSEL",
  var2 = "param", lev2 = "theta",
  var3 = "N_person", lev3 = "N = 100"
)

list_sub



###'#######################################################################
###'
###' `plot_4_factors_list_sub()`
###'
###'

plot_4_factors_list_sub <- function(list_sub, 
                                    y = "mean_log", 
                                    x = "model", 
                                    group = "sum_method", 
                                    facet = "DGM", 
                                    panel = "WLE_rel"){
  
  ### Generate the main plot (without panel)
  p <- ggplot(data = list_sub[[3]], 
              aes(x = .data[[x]], y = .data[[y]], 
                  group = .data[[group]], 
                  color = .data[[group]], 
                  shape = .data[[group]])) +
    
    geom_point(aes(y = .data[[y]]), 
               size = 3, 
               position = position_dodge(width = 0.4)) +  
    
    geom_line(aes(y = .data[[y]]), 
              position = position_dodge(width = 0.4)) +   
    # geom_hline(yintercept = 0, color = "gray30", linetype = "dashed")
    
    geom_text(aes(label = sprintf("%1.2f", round(.data[[y]], 3))),
              hjust = 1.5, color = "black", size = 3,
              position = position_dodge(width = 0.4)) +  
    
    facet_grid(rows = vars(.data[[panel]]),
               cols = vars(.data[[facet]]),
               scales = "free_y") + 
    
    scale_x_discrete(expand = expansion(mult = 0.4)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) + 
    
    # scale_color_manual(values = color_palette[seq(unique(list_sub[[3]]$group))]) +
    # scale_shape_manual(values = shape_palette[seq(unique(list_sub[[3]]$group))]) +
    
    theme_trend +
    theme(strip.background = element_rect(fill = "gray100")) + 
    labs(title = list_sub[[2]])
  
  p
}


### Test
list_sub <- fix_3_factors(
  df_long, 
  var1 = "param", lev1 = "theta",
  var2 = "loss_est", lev2 = "KS_dist",
  var3 = "DGM", lev3 = "Mixed"
)

list_sub

plot_4_factors_list_sub(list_sub, 
                        y = "mean_log",
                        x = "model", 
                        group = "sum_method",
                        facet = "N_person", 
                        panel = "WLE_rel")
