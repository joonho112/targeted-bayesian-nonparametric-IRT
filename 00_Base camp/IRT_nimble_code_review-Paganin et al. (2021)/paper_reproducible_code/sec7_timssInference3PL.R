##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: June, 22 2021
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## nimble version 0.11.1
##-----------------------------------------##
##-----------------------------------------##
# This files compute some quantities for plots and tables
dir.create("figures/dataForFigures", recursive = TRUE, showWarnings = FALSE)
##-----------------------------------------#
## Data
# rm(list = ls()# 
library(bayestestR) ## For HDI intervals
# ## Compute simulation MSE
source("R_functions/ggplot_settings.R")
#############################


##------------------------------------------------------------#
##-----------------------------------------#
library(ggplot2)
library(cowplot)
library(bayestestR)
source("R_functions/ggplot_settings.R")
source("R_functions/multimodalDensity.R")
##-----------------------------------------#
## Set dimensions
plot_width  <- 21*1.1 ## a4 paper
plot_height <- plot_width/5*2
unit   <- "cm" 

##-----------------------------------------#
## Data Timss
##-----------------------------------------#
timssRes <- readRDS("figures/dataForFigures/timss3PL.rds")
##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr <- data.frame(Parametric     = timssRes$paraEstimates$lambda,
                Semiparametric = timssRes$bnpEstimates$lambda)


pTimssDiscr <- ggplot(estimateDiscr, aes(x = Parametric, y =  Semiparametric)) + 
        geom_point(size = 1.5) + 
        scale_x_continuous(breaks = seq(0.5, 3, by = 0.5)) +
        scale_y_continuous(breaks = seq(0.5, 3, by = 0.5)) +
          labs(x = "Parametric", y = "Semiparametric", title = "Discrimination parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pTimssDiscr

ggsave(filename = "figures/3PL_data_timss_discriminations.png", plot = pTimssDiscr,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff <- data.frame(Parametric     = timssRes$paraEstimates$beta,
               Semiparametric = timssRes$bnpEstimates$beta)


pTimssDiff <- ggplot(estimateDiff, aes(x = Parametric, y =  Semiparametric)) + 
        geom_point(size = 1.5) + 
        scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
        scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
          labs(x = "Parametric", y = "Semiparametric", title = "Difficulty parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pTimssDiff

ggsave(filename = "figures/3PL_data_timss_difficulties.png", plot = pTimssDiff,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Guessing
##-----------------------------------------#

estimateGuess <- data.frame(Parametric     = timssRes$paraEstimates$delta,
                           Semiparametric = timssRes$bnpEstimates$delta)


pTimssGuess <- ggplot(estimateGuess, aes(x = Parametric, y =  Semiparametric)) + 
        geom_point(size = 1.5) + 
        scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
        scale_y_continuous(breaks = seq(0, 0.5, by = 0.1)) +
          labs(x = "Parametric", y = "Semiparametric", title = "Guessing parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pTimssGuess

ggsave(filename = "figures/3PL_data_timss_guessing.png", plot = pTimssGuess,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Density - mean 
##-----------------------------------------#

dfEtaMeansTimss <- data.frame(mean = c(timssRes$paraEstimates$eta,
                       timssRes$bnpEstimates$eta))

dfEtaMeansTimss$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansTimss)[1]/2)

title <- paste0("TIMSS data\nDistribution of individual posterior mean abilities")

pEtaMeanTimss <- ggplot(dfEtaMeansTimss, aes(x=mean, color = Model)) +
          geom_histogram(aes(y=..density..), position="identity",
            binwidth = 0.1, fill="white",key_glyph = "path")+
          geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
          xlim(-10, 6) +
            labs(title=title,x="Ability", y = "Density")+
          scale_color_manual(values=c(paraColor, bnpColor),
            guide=guide_legend(override.aes=list(linetype=c(1,1))))


pEtaMeanTimss <- pEtaMeanTimss + theme(legend.title = element_blank())
pEtaMeanTimss

ggsave(filename = "figures/3PL_data_timss_posterior_means.png", plot = pEtaMeanTimss,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid            = timssRes$grid, 
                   semiparametric  = apply(timssRes$densityDPMeasure, 2, mean),
                   parametric      = apply(timssRes$densitySamplesPara, 2, mean))


dfEtaDensityPlot <- reshape2::melt(dfEtaDensity, id.vars = "grid")
colnames(dfEtaDensityPlot) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot$Model) <- c("Semiparametric", "Parametric")
dfEtaDensityPlot$Model <- factor(dfEtaDensityPlot$Model, levels = c("Parametric", "Semiparametric"))

dfIntervalPara <- data.frame(grid  = timssRes$grid, 
                             lower = apply(timssRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(timssRes$densitySamplesPara, 2, quantile, 0.975))
dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = timssRes$grid, 
                             lower = apply(timssRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(timssRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"

title <- paste0("TIMSS data\nEstimate of the distribution of ability")

pEtaDensityTimss <- ggplot(dfEtaDensityPlot, aes(x=grid, y = value, color = Model)) +
            geom_line(lwd = 1, aes(linetype=Model)) +
            geom_line(data = dfIntervalBNP, aes(x = grid, y=lower), lwd = 1,  linetype="dashed") + 
            geom_line(data = dfIntervalBNP, aes(x = grid, y=upper), lwd = 1,  linetype="dashed") +
            geom_line(data = dfIntervalPara, aes(x = grid, y=lower), lwd = 1, linetype="dashed") + 
            geom_line(data = dfIntervalPara, aes(x = grid, y=upper), lwd = 1, linetype="dashed") +
            xlim(-10, 6) +
            scale_linetype_manual(values=c("solid", "solid")) +
            labs(title=title, x="Ability", y = "Density")+
            scale_color_manual(values=c(paraColor, bnpColor),
                guide=guide_legend(override.aes=list(linetype=c("solid", "solid"))))

pEtaDensityTimss <- pEtaDensityTimss + theme(legend.title = element_blank())
pEtaDensityTimss

ggsave(filename = "figures/3PL_data_timss_posterior_density.png", plot = pEtaDensityTimss,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(timssRes$paraPerc, 2, mean),
                                         apply(timssRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(timssRes$paraPerc, 2, quantile, 0.025),
                                         apply(timssRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(timssRes$paraPerc, 2, quantile, 0.975),
                                         apply(timssRes$bnpPerc, 2, quantile, 0.975)))

dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pTimssPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
          width = 0.8, position = position_dodge(width = 0.6)) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(y = "Percentile", x = "Individual", 
                title = "TIMSS data") 

pTimssPerc

ggsave(filename = "figures/3PL_data_timss_percentiles.png", plot = pTimssPerc,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


