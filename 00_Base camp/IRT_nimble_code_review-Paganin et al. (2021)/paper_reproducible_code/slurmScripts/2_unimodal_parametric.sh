#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o unimodal_parametric_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5
######################


FILES=(models/parametric/*.R)

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default


Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered


# Rscript 1_runNimbleModels.R  \
# --model=models/parametric/parametric_IRT_constrainedItem.R \
# --data=data/simulation_unimodal.rds \
# --niter=50000 \
# --nburnin=5000 \
# --nthin=1 \
# --mode=default

# Rscript 1_runNimbleModels.R  \
# --model=models/parametric/parametric_IRT_constrainedAbilities.R \
# --data=data/simulation_unimodal.rds \
# --niter=50000 \
# --nburnin=5000 \
# --nthin=1 \
# --mode=default
