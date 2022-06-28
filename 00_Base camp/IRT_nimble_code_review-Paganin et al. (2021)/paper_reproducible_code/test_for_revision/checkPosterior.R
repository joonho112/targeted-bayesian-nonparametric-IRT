## Check posterior samples 3PNO


res <- readRDS("simulation_bimodal/3PNO/3PNO_mixture.rds")

samples <- res$samples$samples



plot(samples[, "delta[1]"], type = "l")

head(samples)
