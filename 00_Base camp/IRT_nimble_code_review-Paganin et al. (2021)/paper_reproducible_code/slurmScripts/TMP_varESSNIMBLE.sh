#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o essNimble.out            # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
######################

# Rscript TMP_checkMultiESSVariability.R  \
# --niter=20000 \
# --nburnin=2000 

# Rscript TMP_checkMultiESSVariability.R  \
# --niter=40000 \
# --nburnin=4000 

Rscript TMP_checkMultiESSVariability.R  \
--niter=50000 \
--nburnin=5000 

