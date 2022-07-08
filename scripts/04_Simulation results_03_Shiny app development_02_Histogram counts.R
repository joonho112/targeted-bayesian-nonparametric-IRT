
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Simulation Results
###' 
###' Task: Shiny App Development
###'       `Histogram Counts` 
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
###' Load the histogram counts file
###'
###'

load_path <- file.path(data_dir2, "df_histogram_collected_new_true_dens.rds")

df_hist <- read_rds(load_path) 



###'#######################################################################
###'
###' `best_choice_hist()`
###'
###'

best_choice_hist <- function(df_hist = df_hist,
                             param_fix = "theta", 
                             DGM_fix = "Gaussian",
                             N_person_fix = "N = 100",
                             WLE_rel_fix = "WLE rel. = 0.6", 
                             abs_limit = 5){
  
  ### Filter the original data
  df_sub <- df_hist %>%
    filter(sum_method != "true") %>%
    filter(
      param == param_fix, 
      DGM == DGM_fix, 
      N_person == N_person_fix, 
      WLE_rel == WLE_rel_fix
    ) %>%
    filter(
      abs(start) <= abs_limit
    )
  
  # Extract label information
  subtitle_lab <- paste(
    paste0("param: ", param_fix), 
    paste0("true dist.: ", DGM_fix), 
    N_person_fix, 
    WLE_rel_fix,
    collapse = "", sep = ", "
  )
  
  labels <- 
    labs(x = "Estimated latent parameter", 
         y = "Density",
         title = "Scaled empirical distribution function vs. True distribution", 
         subtitle = subtitle_lab, 
         caption = "black solid line: true densities")
  
  ### Plot trellis graph 
  p <- ggplot(df_sub) + 
    geom_rect(aes(xmin = start, xmax = end, 
                  ymin = 0, ymax = density), 
              fill = "limegreen",  color = "gray80", size = 0.0001) + 
    geom_line(aes(x = middle, y = true_dens), size = 0.6) +
    facet_grid(rows = vars(model), cols = vars(sum_method)) + 
    theme_trend + labels
  
  ### Return the resulting object
  list(df_sub, p)
}


### Test the function
best_choice_hist(df = df, 
                 param_fix = "theta", 
                 DGM_fix = "Mixed",
                 N_person_fix = "N = 500",
                 WLE_rel_fix = "WLE rel. = 0.7")


