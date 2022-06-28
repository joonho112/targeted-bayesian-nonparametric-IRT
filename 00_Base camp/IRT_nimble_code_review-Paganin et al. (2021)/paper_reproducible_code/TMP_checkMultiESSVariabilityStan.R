# args <- list()	
# args$nsamples=5000 
# args$nwarmup=5000 

args <- R.utils::commandArgs(asValue=TRUE)

MCMCcontrol <- list()
MCMCcontrol$nwarmup <- as.numeric(args$nwarmup)
MCMCcontrol$niter <- as.numeric(args$nsamples) + MCMCcontrol$nwarmup
##----------------------------------##
##----------------------------------##
## load library and functions
library(rstan)
library(reshape2)
library(mcmcse)
source("R_functions/rescalingFunctions.R")

rstan_options(auto_write = TRUE)

data <- "data/simulation_unimodal_I_10_N_1000.rds"	 

## Load data
if(grepl("timss", data)){
	alldata <- readRDS(data)
	data <- list(y = alldata$y)
} else {
	data <- list(y = readRDS(data))
}

##-----------------------------------------##
## Set variables 
##-----------------------------------------##

##----------------------------------##
cat("Stan 2PL IRT constrained abilities model\n")
##----------------------------------##

fileStan <- "models/parametric_IRT_constrainedAbilities.stan"

## Reshape data in wide format
wide    <- as.data.frame(data)
wide$id <- 1:nrow(wide)  # Attach a person ID number to each row.
long    <- melt(wide, id.vars = "id", variable.name = "item", value.name = "response")
# head(long)

key <- 1:length(unique(long$item))
names(key) <- unique(long$item)
long$item.id <- key[long$item]
# head(long)

stan_data <- list(I  = max(long$item.id), 
				  J  = max(long$id),
				  N  = nrow(long), 
				  ii = long$item.id, 
                  jj = long$id,
                  y  = long$response)

## Create arguments lists
stan_model_args <- list() 
sampling_args   <- list()

## modify stan_model_args
stan_model_args$file <- fileStan
## Create stan_model object
compileTime <- system.time(stan_mod <- do.call(rstan::stan_model, stan_model_args))

## modify sampling args
sampling_args$object <- stan_mod ## object of class stanmodel
sampling_args$data   <- stan_data
sampling_args$chains <- 1

##  Note: in rstan::sampling function the `iter` argument comprises also the number of warmup iterations
sampling_args$warmup <- MCMCcontrol$nwarmup
sampling_args$iter   <- MCMCcontrol$niter 
sampling_args$thin   <- 1
     

##################     
nTimes <- 20
out <- numeric(20)
resList <- list()

paramsList <-list()
adaptList <-list()

for(i in 1:nTimes){
	sampling_args$seed   <- i
	totalTime <- system.time(stan_out <- do.call(rstan::sampling, sampling_args))

	monitors <- c("beta", "lambda", "eta")

	samplesArray <- rstan::extract(stan_out, 
	                            pars = monitors,
	                            permuted = FALSE,
	                            inc_warmup = FALSE)[, 1, ]

 	paramsList[[i]] <- rstan::get_sampler_params(stan_out)
	adaptList[[i]] <- rstan::get_adaptation_info(stan_out)


	modelRes <- posteriorRescalingBeta(samples  = samplesArray[, -grep("^eta", colnames(samplesArray))],
	                                 samples2 = samplesArray[, grep("^eta", colnames(samplesArray))],
	                                 thinEta  = 1, 
	                                 rescale  = TRUE)

	resList[[i]] <- modelRes
	## matrices for ESS evaluations
	onlyItems <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
	              modelRes[grepl("beta", names(modelRes))][[1]])

	itemsAndAbility <- cbind(onlyItems , modelRes$etaSamp)
	itemsAndAbilityMultiESS <- itemsAndAbility[, !grepl("(beta\\[1\\])|(lambda\\[1\\])", colnames(itemsAndAbility))]

	try(out[i] <- mcmcse::multiESS(itemsAndAbilityMultiESS, 
	                                            method = "bm",
	                                            r = 1, 
	                                            adjust = FALSE))
}

saveRDS(out, file = paste0("output/Stan_multiESS_", MCMCcontrol$niter,"_warmup_", MCMCcontrol$nwarmup, ".rds"))

outList <- list(samples = resList, 
				adaptation = adaptList, 
				params = paramsList)

saveRDS(outList, file = paste0("output/Stan_multiESS_res.rds"))

