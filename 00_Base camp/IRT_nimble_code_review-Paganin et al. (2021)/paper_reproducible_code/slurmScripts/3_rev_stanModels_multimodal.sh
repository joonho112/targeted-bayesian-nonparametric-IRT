#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o stan_multimodal.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
###########################


Rscript 1_runStanModel.R  \
--data=data/simulation_multimodal.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_multimodal_I_30_N_1000.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_multimodal_I_30_N_5000.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_multimodal_I_10_N_5000.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_multimodal_I_10_N_1000.rds \
--nsamples=10000 \
--nwarmup=5000 
