##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)

# args <- list()
# args$resFileName <- "output/posterior_samples_submission/simulation_unimodal/parametric/parametric_IRT_unconstrained.rds"


# args$resFileName <- "output/posterior_samples_submission/simulation_unimodal/bnp/bnp_IRT_unconstrained.rds"
# args$elaboratedResFileName <- "output/posterior_samples_elaborated_submission/simulation_unimodal/bnp/bnp_IRT_unconstrained.rds"

## --resFileName
##-----------------------------------------#
data <- strsplit(args$resFileName, "\\/|.rds")[[1]][3]

fileName <- strsplit(args$resFileName, "\\/|.rds")[[1]][5]

## parameterization
param <- strsplit(basename(fileName), "\\_|.rds")[[1]][2]
## model type
modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]
## constraint
constraint <- strsplit(basename(fileName), "\\_|.rds")[[1]][3]

elaboratedResFileName <- gsub("posterior_samples_submission", 
                                "posterior_samples_elaborated_submission",
                                args$resFileName)

## read objects
resObj <- readRDS(args$resFileName)
elaboratedResObj <- readRDS(elaboratedResFileName)
##-----------------------------------------#
## Get efficency results on thinned samples
##-----------------------------------------#
indexThinning <- seq(0, 45000, by = 10)


ess_coda_old   <- NA

## Q: how to consider abilities parameters? 
## rescaling problem
# ess_coda_all   <- NA
ess_coda_latent  <- NA

# ess_multi_all   <- NA
ess_multi_latent  <- NA

onlyItem <- cbind(elaboratedResObj$lambdaSamp, 
            elaboratedResObj$betaSamp[, -1])

latentPars <- cbind(onlyItem[indexThinning, ], 
                    elaboratedResObj$etaSamp)

# latentPars <- cbind(onlyItem, 
#                     elaboratedResObj$etaSamp)

####################
# library(mcmcse)
# multiESS(latentPars[30000:45000, ],  method = "obm")

# cov <- mcse.multi(latentPars, method = "tukey")
# sum(log(cov$eigen_values))

# # multiESS(latentPars, cov = cov, blather = TRUE)


# covmat <- mcse.multi(x, ...)$cov

# var_mat <- cov(chain)
# log.det.var.p <- sum(log(eigen(var_mat, symmetric = TRUE,
# only.values = TRUE)$values))
# log.det.covmat.p <- sum(log(eigs_cov))
# ess <- n * exp((log.det.var.p - log.det.covmat.p)/p)
# return(ess)
####################

## compute ess for item parameters
ess_coda_old <- min(coda::effectiveSize(onlyItem))

ess_coda_latent <- min(coda::effectiveSize(latentPars))
ess_coda_latent_mean <- mean(coda::effectiveSize(latentPars))


# Any of “‘bm’”,“‘obm’”,“‘bartlett’”,“‘tukey’”. “‘bm’”
# represents batch means estimator, “‘obm’” represents the
# overlapping batch means estimator, and “‘bartlett’” and
# “‘tukey’” represent the modified-Bartlett window and the
# Tukey-Hanning windows for the spectral variance estimators.

ess_multi_latent <- mcmcse::multiESS(latentPars)
ess_multi_latent_obm <- mcmcse::multiESS(latentPars, method = "obm")
ess_multi_latent_bartlett <- mcmcse::multiESS(latentPars, method = "bartlett")
ess_multi_latent_tukey <- mcmcse::multiESS(latentPars, method = "tukey")

var_mat <- cov(latentPars)
eigenValsSample <- eigen(var_mat, symmetric = TRUE, only.values = TRUE)
detSampleCov <- sum(log(eigenValsSample$values))



compilationTime <- resObj$compilationTime[3]
runningTime <- resObj$runningTime[3]
samplingTime <- 0

## if parametric save also sampling time
if(modelType == "parametric"){ 
  if(constraint == "stan") { 
    samplingTime <- resObj$samplingTime
  } else {
    percBurnin <- resObj$MCMCcontrol$nburnin/resObj$MCMCcontrol$niter
    samplingTime <- runningTime * (1 - percBurnin)
  }

}

## adding WAIC
WAIC <- 0

if(modelType == "parametric"){ 
  WAIC <- resObj$modelWAIC
} 

outDirTime <- paste0("output/mcmc_time_test/", data)
dir.create(file.path(outDirTime), recursive = TRUE, showWarnings = FALSE)

outFile <- paste0(outDirTime, "/", modelType, "_efficiency.txt")
row <- cbind(fileName, ess_coda_old, ess_coda_latent, ess_coda_latent_mean, ess_multi_latent, ess_multi_latent_obm, 
ess_multi_latent_bartlett, ess_multi_latent_tukey, detSampleCov, compilationTime, runningTime, samplingTime, WAIC)

if(!file.exists(outFile)){
    cat(colnames(row), "\n", file = outFile)
}
# append row
cat(row, "\n", file = outFile, append = TRUE)
# ############################################################################



# res <- readRDS("output/posterior_samples_elaborated_submission/simulation_unimodal/parametric/parametric_IRT_unconstrained.rds")
# res <- readRDS("output/posterior_samples_elaborated_submission/simulation_unimodal/parametric/parametric_IRT_constrainedAbilities.rds")
# res <- readRDS("output/posterior_samples_elaborated_submission/simulation_unimodal/parametric/parametric_IRT_constrainedItem.rds")
# res <- readRDS("output/posterior_samples_elaborated_submission/simulation_unimodal/bnp/bnp_IRT_unconstrained.rds")


# str(res)

# library(mcmcse)
# xx <- cbind(res[grepl("lambda", names(res))][[1]],
#             res[grepl("beta", names(res))][[1]])


# multiESS(xx)

# index <- seq(0, 45000, by = 10)

# allSamp <- cbind(xx[index], res[grepl("^eta", names(res))][[1]])
# multiESS(allSamp)


# cov <- mcse.multi(xx)
# multiESS(x = xx, covmat = cov$cov )

## No need to check on non elaborated samples

# res1 <- readRDS("output/posterior_samples_submission/simulation_unimodal/parametric/parametric_IRT_unconstrained.rds")
# res1 <- readRDS("output/posterior_samples_submission/simulation_unimodal/parametric/parametric_IRT_constrainedAbilities.rds")
# res1 <- readRDS("output/posterior_samples_submission/simulation_unimodal/parametric/parametric_IRT_constrainedItem.rds")

# str(res1)

# index <- seq(0, 45000, by = 10)

# str(res1$samples)

# matSamp <- cbind(res1$samples$samples2, res1$samples$samples[index, ])
# ess_multi <- multiESS(matSamp) 
# ess_multi



#############################
# ## unimodal unconstrained
# ## rescaled
# multiESS(allSamp)
# [1] 7046.4

# ## unimodal IRT constrained abilities
# ## rescaled
# multiESS(allSamp)
# [1] 4500

# ## unimodal IRT constrained items
# ## rescaled
# multiESS(allSamp)
# [1] 7054.657
