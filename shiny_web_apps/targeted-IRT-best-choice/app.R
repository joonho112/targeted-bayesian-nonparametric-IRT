
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Shiny Web App Development
###' 
###' Task: Shiny App Development
###'       `Best Choice` 
###'       
###' Data: 
###' - Performance evaluator estimates (MSEL, MSELP, ISEL, etc)
###' - Histogram Count data
###' 
###' Date: 
###' - 2022-07-07
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###' 

###'#######################################################################
###'
###' Basic settings
###'
###'

### Call libraries
library(shiny)
library(shinythemes)
# library(shinyWidgets)
library(tidyverse)
library(rlang)
library(ggeffects)
library(DT)
library(glue)
library(clubSandwich)


# ### Set temporary working directory for testing
# work_dir <- file.path(path.expand("~"), 
#                       "Documents",
#                       "targeted-bayesian-nonparametric-IRT", 
#                       "shiny_web_apps", 
#                       "targeted-IRT-best-choice") 
# 
# setwd(work_dir)
# 
# list.files(file.path(work_dir, "functions"), full.names = TRUE) %>%
#   walk(source)


###' Call custom functions
list.files(file.path("functions"), full.names = TRUE) %>% 
  walk(source)



###'#######################################################################
###'
###' Prepare datasets
###'
###'

df_long <- reactive({
  
  load_path <- file.path(
    "datasets", 
    "df_loss_estimates_LONG_collected.rds"
  )
  
  read_rds(load_path)
})

df_sum <- reactive({
  load_path <- file.path(
    "datasets", 
    "df_loss_estimates_SUMMARY_collected.rds"
  )
  
  read_rds(load_path)
})

df_hist <- reactive({
  load_path <- file.path(
    "datasets", 
    "df_histogram_collected_new_true_dens.rds"
  )
  
  read_rds(load_path)
})

# df_qt_sum <- reactive({
#   load_path <- file.path(
#     "datasets", 
#     "All_ECDF_estimates_Final_SUMMARY_DP-diffuse_corrected_with_J-300.rds"
#   )
#   
#   read_rds(load_path)
# })  


###'######################################################################
###'
###' User interface
###'
###'

ui <- fluidPage(
  
  ### Shiny theme
  theme = shinytheme("lumen"), 
  
  
  ### Application title
  titlePanel("Targeted Bayesian IRT: Best Model Choice", 
             windowTitle = "Best Model Choice"),
  
  
  ### Sidebar layout with a input and output definitions
  
  sidebarLayout(
    
    ###'############
    ###' Inputs  ###
    ###'############ 
    
    
    sidebarPanel(
      
      ###'#######################################
      ###' Input section 1) Data Generating Conditions
      ###'
      ###' - Select conditions to filter out
      ###'
      
      h3("Data Generating Conditions"),  
      
      # (1) True G distribution
      selectInput(
        inputId = "DGM_fix", 
        label = "True G distribution:", 
        choices = 
          c(
            "Gaussian" = "Gaussian", 
            "Gaussian Mixture (Mixed)" = "Mixed", 
            "Asymmetric Laplace (ALD)" = "ALD"
          ), 
        selected = c("Mixed")
      ),
      
      # (2) N_person: Number of persons 
      selectInput(
        inputId = "N_person_fix", 
        label = "Number of test takers (N):", 
        choices = c("N = 20", 
                    "N = 50", 
                    "N = 100", 
                    "N = 200", 
                    "N = 500"), 
        selected = c("N = 100")
      ),
      
      # (3) WLE reliability
      selectInput(
        inputId = "WLE_rel_fix", 
        label = "WLE reliability:", 
        choices = c("WLE rel. = 0.5", 
                    "WLE rel. = 0.6", 
                    "WLE rel. = 0.7",
                    "WLE rel. = 0.8", 
                    "WLE rel. = 0.9"), 
        selected = c("WLE rel. = 0.7")
      ),

      # # Show data table
      # checkboxInput(inputId = "drop_ML",
      #               label = "Don't show observed ML estimates",
      #               value = TRUE),
      
      hr(),
      
      ###'#######################################
      ###' Input section 2) Choose parameter
      ###'
      
      h3("Parameters"),
      
      selectInput(
        inputId = "param_fix", 
        label = "Parameter:", 
        choices = c("Person Latent Trait (theta)" = "theta",
                    "Item Difficulty (beta)" = "beta"), 
        selected = c("theta")
      ),
      
      
      hr(),
      
      ###'#######################################
      ###' Input section 3) Choose performance evaluator  
      ###'
      
      h3("Performance Evaluator"),
      
      selectInput(
        inputId = "loss_est_fix", 
        label = "Loss function:", 
        # choices = c("MSEL", "MSELR", "MSELP", "ISEL", "IAEL", "KS_dist")
        choices = c("Mean Squared Error Loss (MSEL)" = "MSEL",
                    "MSEL of Rank (MSELR)" = "MSELR",
                    "MSEL of Percentile (MSELP)" = "MSELP",
                    "Integrated Squared Error Loss (ISEL)" = "ISEL",
                    "Integrated Absolute Error Loss (IAEL)" = "IAEL",
                    "Kolmogorov–Smirnov Distance (KS_dist)" = "KS_dist"), 
        selected = c("ISEL")
      ),
      
      
      hr(),
      

      ###'#######################################
      ###' Input section 4) Miscellaneous Inputs
      ###' 
      
      # Built with Shiny by JoonHo Lee
      br(), br(),
      h5("Built with",
         img(src = "https://www.rstudio.com/wp-content/uploads/2014/04/shiny.png", 
             height = "30px"),
         "by JoonHo Lee (jlee296@ua.edu)"
      )
      
      , width = 3),  # End of sidbarPanel
    
    
    
    ###'############
    ###' Outputs
    ###'############ 
    
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  id = "tabsetpanel",
                  
                  tabPanel(title = "Predicted Loss",
                           plotOutput(outputId = "best_choice_cond_pred",
                                      width = "100%")),
                  
                  # tabPanel(title = "Predicted Loss (Table)",
                  #          DT::dataTableOutput(outputId = "df_cond_pred")),
                  
                  tabPanel(title = "Best Choice Table", 
                           br(),
                           DT::dataTableOutput(outputId = "best_choice_table"), 
                           downloadButton(outputId = "download_data", 
                                          label = "Download data")), 
                  
                  tabPanel(title = "Best Choice Map", 
                           plotOutput(outputId = "best_choice_map", 
                                      width = "100%")),
                  
                  tabPanel(title = "Scaled EDF", 
                           plotOutput(outputId = "best_choice_hist", 
                                      width = "100%"))
                  
                  # tabPanel(title = "Bias in Percentile", 
                  #          plotOutput(outputId = "best_choice_percentile", 
                  #                     width = "100%"))
                  
      ) # End of tabsetPanel
      , width = 9)  # End of mainPanel
    
  )  # End of sidebarLayout
)  # End of ui() function



###'######################################################################
###'
###' Server logic
###'
###'

server <- function(input, output, session) {
  
  ###'###################################################
  ###' Create reactive data table - conditional predictions
  ###' 
  
  df_cond_pred <- reactive({
    
    # Filter/Subset the data
    data <- df_long() %>%
      filter(
        loss_est == input$loss_est_fix, 
        param == input$param_fix, 
        DGM == input$DGM_fix, 
        N_person == input$N_person_fix, 
        WLE_rel == input$WLE_rel_fix
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
      mutate(loss_est = input$loss_est_fix) %>%
      dplyr::select(loss_est, everything())
    
    df_cond_pred
  })
  
  ###'###################################################
  ###' Create reactive data table: Predicted Loss (Table)
  ###' 
  
  tbl_cond_pred <- reactive({
    
    df_cond_pred() %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      dplyr::select(loss_est, x, group, 
                    predicted, std.error, conf.low, conf.high) %>%
      arrange(predicted) %>%
      set_names("loss", "model", "summary",  
                "model-predicted loss", "standard error", 
                "lower_95CI", "upper_95CI")
  })
  
  
  # output$df_cond_pred <- DT::renderDataTable({
  #     
  #     tbl_cond_pred <- df_cond_pred() %>%
  #         mutate(across(where(is.numeric), ~round(., 3))) %>%
  #         dplyr::select(loss_est, x, group, 
  #                       predicted, std.error, conf.low, conf.high) %>%
  #         arrange(predicted) %>%
  #         set_names("loss", "model", "summary",  
  #                   "model-predicted loss", "standard error", 
  #                   "lower_95CI", "upper_95CI")
  #         
  #     DT::datatable(data = tbl_cond_pred, 
  #                   options = list(pageLength = 12), 
  #                   rownames = FALSE)
  # })
  
  
  ###'###################################################
  ###' Create Output: Predicted Loss (Plot)
  ###' 
  
  output$best_choice_cond_pred <- renderPlot({
    
    # Extract label information
    subtitle_lab <- paste(
      paste0("parameter: ", input$param_fix), 
      paste0("True G: ", input$DGM_fix), 
      input$N_person_fix, 
      input$WLE_rel_fix,
      collapse = "", sep = ", "
    )
    
    labels <- 
      labs(x = "Model for G", 
           y = paste0("Predicted (", input$loss_est_fix, ")"), 
           title = paste0("Predicted ",  input$loss_est_fix," by meta−model regression"),
           subtitle = subtitle_lab, 
           caption = NULL)
    
    # Create a plot
    # Generate the main plot (without panel)
    p <- ggplot(data = df_cond_pred(),
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
      
      scale_color_manual(values = color_palette[seq(unique(df_cond_pred()$group))]) +
      scale_shape_manual(values = shape_palette[seq(unique(df_cond_pred()$group))]) +  
      
      theme_trend +
      theme(strip.background = element_rect(fill = "gray100")) + 
      labels
    
    # Return only the plot
    p
    
  }, width = 800, height = 600)
  
  
  
  ###'###################################################
  ###' Create Output: Best Choice Table
  ###' 
  
  ### Define reactive data table
  tbl_best_choice <- reactive({
    
    df_temp <- best_choice_table(df_sum = df_sum(),
                                 param_fix = input$param_fix, 
                                 DGM_fix = input$DGM_fix,
                                 N_person_fix = input$N_person_fix,
                                 WLE_rel_fix = input$WLE_rel_fix) %>%
      filter(loss_est == input$loss_est_fix) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      dplyr::select(loss_est, model, sum_method, 
                    mean_raw, delta_raw, pct_delta_raw) %>%
      set_names("loss", "model", "summary",
                "mean loss", "difference", "percentage change") %>%
      right_join(tbl_cond_pred(), by = c("loss", "model", "summary")) %>%
      mutate(rank = 1:n()) %>%
      relocate(rank, .after = "loss")
    
    df_temp
  })

  ### Print data table
  output$best_choice_table <- DT::renderDataTable(
    
    DT::datatable(data = tbl_best_choice(), 
                  options = list(pageLength = 12), 
                  rownames = FALSE)
  )
  
  ### Download file
  output$download_data <- downloadHandler(
    
    filename = function() {
      
      paste0("data_table.csv")
      
    },
    
    content = function(file) { 
      
      write.csv(tbl_best_choice(), file) 
      
    }
  )
  
  
  
  ###'###################################################
  ###' Create Output: (3) Scaled EDF
  ###' 
  
  output$best_choice_hist <- renderPlot({
    
    best_choice_hist(df = df_hist(), 
                     param_fix = input$param_fix, 
                     DGM_fix = input$DGM_fix,
                     N_person_fix = input$N_person_fix,
                     WLE_rel_fix = input$WLE_rel_fix, 
                     abs_limit = 5)
    
  }, width = 1000, height = 600)
  
  
  ###'###################################################
  ###' Create Output: (4) Best choice heatmap
  ###' 
  
  output$best_choice_map <- renderPlot({
    
    best_choice_map(df = df_sum(), 
                    loss_est_fix = input$loss_est_fix, 
                    param_fix = input$param_fix, 
                    DGM_fix = input$DGM_fix,
                    N_person_fix = input$N_person_fix,
                    WLE_rel_fix = input$WLE_rel_fix)
    
  }, width = 800, height = 600)
  
  
  # ###'###################################################
  # ###' Create Output: (5) Biases in percentile estimates
  # ###' 
  # 
  # output$best_choice_percentile <- renderPlot({
  #   
  #   best_choice_percentile(df = df_qt_sum(), 
  #                          true_G_fix = input$true_G_fix,
  #                          J_fix = input$J_fix,
  #                          sigma_tau_fix = input$sigma_tau_fix,
  #                          nj_mean_fix = input$nj_mean_fix,
  #                          cv_fix = input$cv_fix, 
  #                          drop_ML = input$drop_ML)
  #   
  # }, width = 1000, height = 650)
  
} # End of server() function



###'######################################################################
###' 
###' Run the application
###' 
###'   
# Run the application 
shinyApp(ui = ui, server = server)

