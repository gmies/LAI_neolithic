#!/bin/bash

#1/9/26 update to run LAI on X 


#11/12/25 update to run on differing sample sizes of sources and compare GA correlations (and LA as well)

#script to run all lai methods on a dataset as an input 
#start this script in the results section of a directory


#script 6/11/25 to rerun with gens 35 for all methods 


#things you add to the script and can edit when you run, since these will change

#put input of data here as an argument to be used in the script
vcf_dir="/project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/data"
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

num_hg = 48
#num_hg=7
num_neo=7


mkdir input_files
cp /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/* input_files/.


#make directory for ancestry_hmm
mkdir ancestry_hmm
cd ancestry_hmm


#make input files:

#make frequency for each ancestry:
module load plink/1.90Beta6.18

#admixed / mneo file:
KEEP_FILE="../input_files/mneo_keep_fam.txt"  # Path to the mneo_keep_fam.txt file

# Initialize output files for the concatenated results
output_neo="neo_combined.frq"
output_hg="hg_combined.frq"
output_bim="combined.bim"

# Clear the output files if they exist
> $output_neo
> $output_hg
> $output_bim

chr=23
# Process each chromosome
    # Run plink to generate bed file
    plink --vcf $vcf_dir/$vcf_name${chr}.vcf.gz --const-fid --make-bed --out full_dataset_chr${chr}

#edit fam file to have sex:

cp full_dataset_chr23.fam saved_full_dataset_chr23.fam
awk 'NR==FNR {sex[$1]=$2; next} {if($2 in sex) $5=sex[$2]; print}' /project/mathilab/gmies/neolithic_selection/X_chromosome/imputation_2025/x_chr_bams/names_ploidy.txt full_dataset_chr23.fam > full_dataset_chr23_updated.fam
mv full_dataset_chr23_updated.fam full_dataset_chr23.fam

    
    # Generate frequency counts for neo dataset
    plink --bfile full_dataset_chr${chr} --keep ../input_files/neo_keep_fam.txt --freq counts --out ${chr}_neo.frq
    
    # Generate frequency counts for hg dataset
    plink --bfile full_dataset_chr${chr} --keep ../input_files/hg_keep_fam.txt --freq counts --out ${chr}_hg.frq
    
    # Concatenate neo frequency files, keeping the header only from the first file
    if [ $chr -eq 1 ]; then
        cat ${chr}_neo.frq.frq > $output_neo
    else
        tail -n +2 ${chr}_neo.frq.frq >> $output_neo
    fi
    
    # Concatenate hg frequency files, keeping the header only from the first file
    if [ $chr -eq 1 ]; then
        cat ${chr}_hg.frq.frq > $output_hg
    else
        tail -n +2 ${chr}_hg.frq.frq >> $output_hg
    fi

    # Append the bim file
    cat full_dataset_chr${chr}.bim >> $output_bim
    
    
    
    #now do for each admixed indiv (mneo)
    # Counter for sequential numbering
    count=1

    # Loop through each line in mneo_keep_fam.txt
    while IFS= read -r line; do
        # Create a temporary keep file for this iteration
        keep_file="${count}_mneo_keep_fam.txt"
        echo "$line" > "$keep_file"

        # Run PLINK command to calculate allele frequencies for the current chromosome
        plink --bfile "full_dataset_chr${chr}" --keep "$keep_file" --freq counts --out "${chr}_${count}_mneo"

        # Remove temporary keep file
        rm "$keep_file"

        # Increment the counter
        ((count++))
    done < "$KEEP_FILE"
    
    

#concat files together:
module load R/4.4

Rscript ../../scripts/ancestry_hmm/concat.R


#combine individuals for each file:

Rscript ../../scripts/ancestry_hmm/combine_indivs.R $num_admixed

#going to delete intermediate files here, so delete the files made from run_chr_mneo.sh
rm *.nosex 
rm *.log
rm *_mneo.frq.counts 



#add cm distance to bim file:

Rscript ../../scripts/ancestry_hmm/X_add_bims_cm.R ./




#know hg is the minor ancestry here (keeping the files the same)


#make frequency for each ancestry:

#make input file with cm distance and no NAs

Rscript ../../scripts/ancestry_hmm/make_input_file.R ./hg_concat.frq ./neo_concat.frq ./ mneo.counts final_allentoft_bim.bim $num_admixed output_hmm_run_file

#edit samples file:
awk 'NR==FNR {sex[$1]=$2; next} {if($1 in sex) $2=sex[$1]; print $1, $2}' /project/mathilab/gmies/neolithic_selection/X_chromosome/imputation_2025/x_chr_bams/names_ploidy.txt ../input_files/ancestry_hmm_samples.txt > ../input_files/ancestry_hmm_samples_ploidy.txt



#run method:

#run ancestry hmm:
module load Ancestry_HMM/1.0.2


ancestry_hmm -i output_hmm_run_file -s ../input_files/ancestry_hmm_samples_ploidy.txt -a 2 0.2 0.8 -p 0 0 0.2 -p 1 $input_gen 0.8 -g --ne 10000




#make output file:
mkdir output
cd output

Rscript ../../../scripts/ancestry_hmm/ploidy_output_average_ancestry_hmm.R ../ output_ancestry_hmm.txt

cp ../output_ancestry_hmm.txt .


#run sliding bins z scores:
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 ../output_ancestry_hmm.txt

cd ..
cd ..
