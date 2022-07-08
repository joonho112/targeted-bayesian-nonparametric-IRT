
###'######################################################################
###'
###' Project: Targeted Bayesian nonparametric IRT
###' 
###' Category: Shiny Web App Development
###' 
###' Task: Define general helper functions
###'       
###' Data: No data. 
###' 
###' Date: 
###' - 2022-07-07
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###' 

###'######################################################################
###'
###' Common settings for theme and temporary labels
###'
###'

### Theme settings
theme_trend <- 
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())

### Temporary labels
temp_labels <- labs(title = "Enter title here", 
                    subtitle = "Enter subtitle here", 
                    caption = "Enter caption here", 
                    y = "Enter ylabel here",  
                    x = "Enter xlabel here")

### Define manual palettes
color_palette <- c("firebrick1", "dodgerblue1", "forestgreen", "darkorchid1", "darkgoldenrod1", 
                   "blue", "green", "purple", "gold", "red")    
shape_palette <- c(16, 17, 15, 18, 1, 2, 0, 5, 6, 4, 3, 8, 10, 7, 9) 


