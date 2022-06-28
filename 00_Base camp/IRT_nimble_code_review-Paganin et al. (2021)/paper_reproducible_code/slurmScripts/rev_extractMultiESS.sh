## Extract results for parametric models 

for filename in output/posterior_samples_submission/simulation_unimodal/parametric/*.rds; do
	echo $filename
	Rscript rev_testmultiESS.R --resFileName=$filename
done

for filename in output/posterior_samples_submission/simulation_bimodal/parametric/*.rds; do
	echo $filename
	Rscript rev_testmultiESS.R --resFileName=$filename
done


for filename in output/posterior_samples_submission/simulation_unimodal/bnp/*.rds; do
	echo $filename
	Rscript rev_testmultiESS.R --resFileName=$filename
done


for filename in output/posterior_samples_submission/simulation_bimodal/bnp/*.rds; do
	echo $filename
	Rscript rev_testmultiESS.R --resFileName=$filename
done
