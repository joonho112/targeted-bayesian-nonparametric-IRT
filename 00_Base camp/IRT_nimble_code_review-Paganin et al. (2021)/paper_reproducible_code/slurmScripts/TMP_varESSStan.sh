#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o essStan3.out            # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
######################

#  Rscript TMP_checkMultiESSVariabilityStan.R  \
# --nsamples=5000 \
# --nwarmup=5000


#  Rscript TMP_checkMultiESSVariabilityStan.R  \
# --nsamples=10000 \
# --nwarmup=10000

#  Rscript TMP_checkMultiESSVariabilityStan.R  \
# --nsamples=50000 \
# --nwarmup=5000


 Rscript TMP_checkMultiESSVariabilityStan.R  \
--nsamples=10000 \
--nwarmup=5000
