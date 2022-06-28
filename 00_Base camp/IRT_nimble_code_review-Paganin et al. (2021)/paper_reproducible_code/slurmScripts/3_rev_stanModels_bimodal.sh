#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o stan_bimodal.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
###########################

Rscript 1_runStanModel.R  \
--data=data/simulation_bimodal.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_bimodal_I_30_N_1000.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_bimodal_I_30_N_5000.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_bimodal_I_10_N_5000.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_bimodal_I_10_N_1000.rds \
--nsamples=10000 \
--nwarmup=5000 
