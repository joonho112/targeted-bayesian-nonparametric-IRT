##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: March 2022
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## nimble version 0.11.1
##-----------------------------------------#
source("R_functions/ggplot_settings.R")
library(cowplot)
##-----------------------------------------#
## Set dimensions
plot_width  <- 21 ## a4 paper
plot_height <- plot_width/5*2
unit   <- "cm" 

dir <- "output/posterior_samples/mcmc_time/"
paraFileName <- "parametric_efficiency.txt"
bnpFileName <- "bnp_efficiency.txt"
##-----------------------------------------#
## Plot efficiency results for simulation
##-----------------------------------------#
fileList <- list.files("output/mcmc_time", full.names = TRUE)

##-----------------------------------------#
## Unimodal - extraSimulation
##-----------------------------------------#

unimodalFiles <- fileList[grep("unimodal", fileList)]

unimodalList <- list()
for(i in 1:length(unimodalFiles)){
      unimodalList[[i]] <- read.table(paste0(unimodalFiles[i], "/", paraFileName), header = TRUE )
      unimodalList[[i]]$simulation <- strsplit(unimodalFiles[i], "\\/")[[1]][3]
}     


unimodalDf <- as.data.frame(do.call(rbind, unimodalList))

# unique(unimodalDf$fileName)
# unique(labelData$R_label)

# unimodalDf$ESS_second <- unimodalDf$essCodaLogLik/unimodalDf$runningTime
# unimodalDf$ESS_second <- unimodalDf$essCodaLogPostItemsAbility/unimodalDf$runningTime
unimodalDf$ESS_second <- unimodalDf$multiEssItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("parametric_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulatio, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("mESS/second (total time)")

pUni <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,110)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Unimodal Scenario")

pUni

ggsave(filename = "figures/SM_unimodalMultiESS.png", plot = pUni,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## Bimodal
##-----------------------------------------#

bimodalFiles <- fileList[grep("bimodal", fileList)]

bimodalList <- list()
for(i in 1:length(bimodalFiles)){
      bimodalList[[i]] <- read.table(paste0(bimodalFiles[i], "/", paraFileName), header = TRUE )
      bimodalList[[i]]$simulation <- strsplit(bimodalFiles[i], "\\/")[[1]][3]

}     

bimodalDf <- as.data.frame(do.call(rbind, bimodalList))

bimodalDf$ESS_second <- bimodalDf$multiEssItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("parametric_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("mESS/second (total time)")

pBi <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,110)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Bimodal Scenario")

pBi

ggsave(filename = "figures/SM_bimodalMultiESS.png", plot = pBi,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", paraFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$multiEssItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("parametric_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])


## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])


pMulti <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,110)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Multimodal Scenario")

pMulti

ggsave(filename = "figures/SM_multimodalMultiESS.png", plot = pMulti,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')
##-----------------------------------------#


allPlots <- plot_grid(
    plot_grid(
    pUni  + theme(legend.position = "none"),
    pBi  + theme(legend.position = "none"),
    pMulti  + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pUni + theme(legend.position = "bottom")), 
   rel_heights = c(1, .1), nrow=2)

allPlots

ggsave(filename = "figures/SM_fig1_allScenarioMultiESS.png", 
        plot = allPlots,
        width = plot_width, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')


##-----------------------------------------#
##-----------------------------------------#
## min ESS 
##-----------------------------------------#

unimodalFiles <- fileList[grep("unimodal", fileList)]

unimodalList <- list()
for(i in 1:length(unimodalFiles)){
      unimodalList[[i]] <- read.table(paste0(unimodalFiles[i], "/", paraFileName), header = TRUE )
      unimodalList[[i]]$simulation <- strsplit(unimodalFiles[i], "\\/")[[1]][3]
}     


unimodalDf <- as.data.frame(do.call(rbind, unimodalList))

unimodalDf$ESS_second <- unimodalDf$essCodaItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("parametric_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulatio, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("minESS/second (total time)")

dfParametricEff$Simulation

pUni <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,5)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Unimodal Scenario")

pUni

ggsave(filename = "figures/SM_unimodalminESS.png", plot = pUni,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## Bimodal
##-----------------------------------------#

bimodalFiles <- fileList[grep("bimodal", fileList)]

bimodalList <- list()
for(i in 1:length(bimodalFiles)){
      bimodalList[[i]] <- read.table(paste0(bimodalFiles[i], "/", paraFileName), header = TRUE )
      bimodalList[[i]]$simulation <- strsplit(bimodalFiles[i], "\\/")[[1]][3]

}     

bimodalDf <- as.data.frame(do.call(rbind, bimodalList))

bimodalDf$ESS_second <- bimodalDf$essCodaItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("parametric_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("minESS/second (total time)")

pBi <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,5)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Bimodal Scenario")

pBi

ggsave(filename = "figures/SM_bimodalminESS.png", plot = pBi,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')


##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", paraFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$essCodaItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("parametric_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])


## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])


colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)
## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

pMulti <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,5)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Multimodal Scenario")

pMulti

ggsave(filename = "figures/SM_multimodalminESS.png", plot = pMulti,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')
##-----------------------------------------#


allPlots <- plot_grid(
    plot_grid(
    pUni  + theme(legend.position = "none"),
    pBi  + theme(legend.position = "none"),
    pMulti  + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pUni + theme(legend.position = "bottom")), 
   rel_heights = c(1, .1), nrow=2)

allPlots

ggsave(filename = "figures/SM_allScenarioMinESS.png", 
        plot = allPlots,
        width = plot_width, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')


