#!/bin/bash

#1/9/26 update to run LAI on X 


#11/12/25 update to run on differing sample sizes of sources and compare GA correlations (and LA as well)

#script to run all lai methods on a dataset as an input 
#start this script in the results section of a directory


#script 6/11/25 to rerun with gens 35 for all methods 


#things you add to the script and can edit when you run, since these will change

#put input of data here as an argument to be used in the script
vcf_dir="/project/mathilab/gmies/neolithic_selection/X_chromosome/imputation_2025/imputation/step3_phased"
#vcf names
vcf_name="chrX.phased."
#and within this dir should be vcfs that look like this: $vcf_name${chr}.vcf.gz

# Set input_gen to the value you want, this is used in: ancestry_hmm
input_gen=35 

#set number of admixed individuals, used in: ancestry_hmm
num_admixed=176

#give name of sample file for admixed samples, in data/samples directory outside of results
admixed_samples=mneo_samples


#saved num ssa which is $num_admixed * 2
num_admixed_count=$((num_admixed * 2))

#set number of admixed individuals, used in: simplai


#hg is ceu and neo is yri


#things you need to manually make and have in the input_files directory in results:

#1. for ancestry hmm: the sample file which is a list of samples (admixed indivs) and then a col of 2 or ploidy #
#sample file example: awk '{print $1, "2"}' /project/mathilab/gmies/neolithic_selection/1kg_runs/data/samples/asw_samples > ancestry_hmm_samples.txt

#2. for ancestry hmm: for minor ancestry: neo_keep_fam.txt and for major ancestry: hg_keep_fam.txt which are 0 and then the sample names, and do the same for admixed as mneo_keep_fam.txt
#examples: awk '{print "0", $1}' /project/mathilab/gmies/neolithic_selection/1kg_runs/data/samples/7_yri_samples > neo_keep_fam.txt
		# awk '{print "0", $1}' /project/mathilab/gmies/neolithic_selection/1kg_runs/data/samples/48_ceu_samples > hg_keep_fam.txt
		# awk '{print "0", $1}' /project/mathilab/gmies/neolithic_selection/1kg_runs/data/samples/asw_samples > mneo_keep_fam.txt

#3. for flare, make reference panel which is 7 neo and neo and then 48 hg (minor ancestry) and hg
#example: awk '{print $1, "neo"}' /project/mathilab/gmies/neolithic_selection/1kg_runs/data/samples/7_yri_samples > flare.ref.panel
		# awk '{print $1, "hg"}' /project/mathilab/gmies/neolithic_selection/1kg_runs/data/samples/48_ceu_samples >> flare.ref.panel


#all scripts ran within will be here: /project/mathilab/gmies/neolithic_selection/1kg_runs/scripts/ and within their method directory


#11/12/25 make input files:

#/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes

#want options of 7v7, 7v3, 7v1
#first number is farmer, second is hg, originally was 7 v 48

num_hg=48
#num_hg=7
num_neo=7


#run mosaic:
#run method:

cd mosaic_7v7

module load R/4.4

awk 'BEGIN{OFS="\t"} {
    gsub(/^X_/, "23_", $1)
    $2 = 23
    print
}' snpfile.23 > snpfile.23.fixed

mv snpfile.23.fixed snpfile.23



awk '{
  for (i=1; i<=NF; i++) {
    if ($i=="NA") $i=0
  }
  print
}' rates.23 > rates.23.fixed

mv rates.23.fixed rates.23



head -160149 snpfile_full.23 > snpfile_test.23
head -160149 rates_full.23 > rates_test.23
head -160149 mneogenofile_full.23 > mneogenofile_test.23
head -160149 neogenofile_full.23 > neogenofile_test.23
head -160149 hggenofile_full.23 > hggenofile_test.23

# Rename original files
mv snpfile.23 snpfile_full.23
mv rates.23 rates_full.23
mv mneogenofile.23 mneogenofile_full.23
mv neogenofile.23 neogenofile_full.23
mv hggenofile.23 hggenofile_full.23

# Rename test files
mv snpfile_test.23 snpfile.23
mv rates_test.23 rates.23
mv mneogenofile_test.23 mneogenofile.23
mv neogenofile_test.23 neogenofile.23
mv hggenofile_test.23 hggenofile.23


Rscript /project/mathilab/gmies/MOSAIC-master/mosaic.R mneo ./ -c 23:23 -n $num_admixed -N 10000 --gens $input_gen 


#make output file:


Rscript ../scripts/mosaic/mosaic_output_file.R ./MOSAIC_RESULTS/localanc_mneo_2way_*.RData ./MOSAIC_RESULTS/mneo_2way_*.RData output_mosaic.txt  















