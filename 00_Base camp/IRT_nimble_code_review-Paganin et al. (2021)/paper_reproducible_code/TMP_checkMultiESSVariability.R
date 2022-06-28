####################################
####################################
## Run a model multiple times and saves the multiESS
## using the small simulation

library(nimble)
library(mcmcse)
source("R_functions/rescalingFunctions.R")

args <- R.utils::commandArgs(asValue=TRUE)


MCMCcontrol 		<- list()
MCMCcontrol$niter 	<- as.numeric(args$niter)
MCMCcontrol$nburnin <- as.numeric(args$nburnin)

# ## MCMC settings
# MCMCcontrol 		<- list()
# MCMCcontrol$niter 	<- 10000
# MCMCcontrol$nburnin <- 1000


data 	<- list(y = readRDS("data/simulation_unimodal_I_10_N_1000.rds"))

source("models/parametric/parametric_IRT_constrainedAbilities.R")

scores 		<- apply(data$y, 1, sum)
Sscores 	<- (scores - mean(scores))/sd(scores)
inits$eta 	<- Sscores


model <- nimbleModel(code = code,
					 data 		= data,  
					 constants	= constants,
					 inits 		= inits, 
					 calculate 	= FALSE)


## update monitors
mcmcConf <- configureMCMC(model, monitors = monitors)
mcmcConf$addMonitors("eta")

mcmc <- buildMCMC(mcmcConf)	




compilationTime <- system.time({
    Cmodel <- try(compileNimble(model))
    if(inherits(Cmodel, 'try-error')) {
      stop("There was a problem compiling the nimble model.")
    }
    Cmcmc <- try(compileNimble(mcmc, project = model))
    if(inherits(Cmodel, 'try-error')) {
      stop("There was a problem compiling the nimble MCMC.")
    }
})


##---------------------------------------------------##
## Run MCMC 
##---------------------------------------------------##
nTimes <- 20
out <- numeric(20)
resList <- list()

for(i in 1:nTimes){

	runningTime <- system.time({try(
		res <- runMCMC(Cmcmc, 
					   niter 	= MCMCcontrol$niter, 
					   nburnin  = MCMCcontrol$nburnin, 
					   setSeed = i))
		if(inherits(res, 'try-error')) {
	  		warning(paste0("There was a problem running nimble MCMC."))
	  }
	})


	modelRes <- posteriorRescalingBeta(samples  = res[, -grep("^eta", colnames(res))],
	                                   samples2 = res[, grep("^eta", colnames(res))],
	                                   thinEta  = 1, 
					                  rescale  = TRUE)

	resList[[i]] <- modelRes

	## matrices for ESS evaluations
	onlyItems <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
	              modelRes[grepl("beta", names(modelRes))][[1]])

	itemsAndAbility <- cbind(onlyItems, modelRes$etaSamp)
	itemsAndAbilityMultiESS <- itemsAndAbility[, !grepl("(beta\\[1\\])|(lambda\\[1\\])", colnames(itemsAndAbility))]

	try(out[i] <- mcmcse::multiESS(itemsAndAbilityMultiESS, 
	                                            method = "bm", r = 1, 
	                                            adjust = FALSE))
}

save(out, file = paste0("output/multiESS_", MCMCcontrol$niter, ".rds"))
saveRDS(resList, file = paste0("output/NIMBLE_multiESS", MCMCcontrol$niter, "res.rds"))



