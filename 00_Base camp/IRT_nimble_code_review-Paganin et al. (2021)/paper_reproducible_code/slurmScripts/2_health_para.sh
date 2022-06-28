#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o health_parametric_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5
######################

FILES=(models/parametric/*.R)

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--dirResults=/scratch/users/sallypaganin \
--data=data/data_health.rds \
--niter=50000 \
--nburnin=25000 \
--nthin=1 \
--mode=default

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--dirResults=/scratch/users/sallypaganin \
--data=data/data_health.rds \
--niter=50000 \
--nburnin=25000 \
--nthin=1 \
--mode=centered