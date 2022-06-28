#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o stan_bimodal.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
###########################

Rscript 1_runStanModel.R  \
--data=data/data_timss.rds \
--nsamples=15000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/data_health.rds \
--nsamples=15000 \
--nwarmup=5000 