#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o multimodal_extra_bnp_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3
######################

FILES=(models/bnp/*.R)

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_multimodal2_I_10_N_1000.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_multimodal2_I_10_N_1000.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_multimodal2_I_10_N_5000.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_multimodal2_I_10_N_5000.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered
