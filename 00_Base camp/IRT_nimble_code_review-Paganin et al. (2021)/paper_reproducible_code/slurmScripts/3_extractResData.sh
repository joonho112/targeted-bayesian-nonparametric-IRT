#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o extractParaBNP3PL.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1

for filename in /scratch/users/sallypaganin/data_timss/parametric3PL/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename 
done
##################################################################
## This has run
# for filename in /scratch/users/sallypaganin/data_health/parametric/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename
# done

# # DONE
# for filename in /scratch/users/sallypaganin/data_health/bnp/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename
# done

############################################################
## Extract results for semiparametric models 

# for filename in /scratch/users/sallypaganin/data_timss/bnp/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename
# done

# Extract results for parametric models 
# for filename in /scratch/users/sallypaganin/data_timss/parametric/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename 
# done

################

# rsync -chavzP --stats sallypaganin@arwen.berkeley.edu:/scratch/users/sallypaganin/