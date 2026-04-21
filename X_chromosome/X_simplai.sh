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
#num_admixed_count=$((num_admixed * 2))

#set number of admixed individuals, used in: simplai

chr=23

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

#num_hg=48
num_hg=7
num_neo=7


#mkdir input_files
#cp /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/* input_files/.
#make directory for ancestry_hmm

#make frequency for each ancestry:
module load plink/1.90Beta6.18


#concat files together:
module load R/4.4


module load bcftools/1.20 

#need number to be females *2 + males

#num_admixed_count=257
num_admixed_count=352

#and run simplai 7v7
#i think just edit to grab the first 7 of the hg (ceu)
#make directory for simplai
mkdir simplai
cd simplai

#make input files:

#this is edited by number of individuals, would need to change for more individuals ***
#and need num snps in mosaic output ***
num_snps=$(wc -l ../mosaic_7v7/mneogenofile.23 | grep -v total | awk '{s+=$1} END {print s}')

Rscript ../scripts/simplai/make_input_simpLAI.R ../mosaic_7v7/ $num_snps $num_hg $num_neo $num_admixed input_simplai.gen

#run method:

#saved num ssa which is $num_admixed * 2

# Ensure script stops on error
set -e

# Define the genome lengths for each chromosome (you can modify these based on actual genome sizes)
declare -A genome_lengths=(
    [23]=155270560
)

i=23
    # Create output directory for each chromosome
    mkdir -p ${i}_chr_output
    cd ${i}_chr_output

   # Extract and write the header (first line) to the output file
    head -n 1 ../input_simplai.gen > chr_${i}.gen
    
    # Extract data for the specific chromosome
    awk -v i="$i" '$1 == i' ../input_simplai.gen >> chr_${i}.gen
    
    # Run simpLAI with the extracted data, using the chromosome-specific genome length
    #/project/mathilab/gmies/neolithic_selection/simpLAI/simpLAI -g chr_${i}.gen --ss1 70 --ss2 11 --ssa $num_admixed_count -l ${genome_lengths[$i]} -s 1e6 -i 5e5 -n 2000 -m 1000 -t 5
    /project/mathilab/gmies/neolithic_selection/simpLAI/simpLAI -g chr_${i}.gen --ss1 14 --ss2 14 --ssa $num_admixed_count -l ${genome_lengths[$i]} -s 1e6 -i 5e5 -n 2000 -m 1000 -t 5
    # Go back to the parent directory
    cd ..



#make output file:



head -n 1 ../ancestry_hmm/output/output_ancestry_hmm.txt > full_output_simpLAI.txt

Rscript ../scripts/simplai/format_simpLAI_output.R ${i}_chr_output/*chr_${i}*fromMin_withSingl.adm $num_admixed_count

mv simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
awk -v i="$i" 'NR == 1 {print $0} NR > 1 {$1 = i; print}' ${i}_simpLAI_final_output.txt > temp_${i}_simpLAI_final_output.txt
mv temp_${i}_simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
awk -v i="$i" '$1 == i' ${i}_simpLAI_final_output.txt >> full_output_simpLAI.txt


cp full_output_simpLAI.txt output/output_simplai.txt

mkdir output
cd output
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 output_simplai.txt

cd ..













