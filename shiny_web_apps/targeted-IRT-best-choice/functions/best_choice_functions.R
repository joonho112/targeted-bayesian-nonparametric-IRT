###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Shiny Web App Development
###' 
###' Task: Define functions for the best choice app
###'       
###' Data: No data. 
###' 
###' Date: 
###' - 2022-07-07
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###' 

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

# ### Test the function
# best_choice_map(df_sum = df_sum, 
#                 param_fix = "theta", 
#                 DGM_fix = "Mixed",
#                 N_person_fix = "N = 100",
#                 WLE_rel_fix = "WLE rel. = 0.9", 
#                 loss_est_fix = "IAEL")



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


# ### Test
# df_best <- best_choice_table(
#   df_sum = df_sum,
#   # loss_est_fix = "ISEL", 
#   param_fix = "theta", 
#   DGM_fix = "Mixed",
#   N_person_fix = "N = 100",
#   WLE_rel_fix = "WLE rel. = 0.6")
# 
# df_temp <- df_best %>%
#   mutate(across(where(is.numeric), ~round(., 3)))
# 
# View(df_temp)



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

# ### Test
# best_choice_cond_pred(df_long = df_long, 
#                       loss_est_fix = "KS_dist", 
#                       param_fix = "theta", 
#                       DGM_fix = "Mixed",
#                       N_person_fix = "N = 200",
#                       WLE_rel_fix = "WLE rel. = 0.7")



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


