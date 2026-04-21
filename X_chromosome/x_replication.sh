
#!/bin/bash
#script to run all lai methods on a dataset as an input 
#start this script in the results section of a directory


#r scripts are in run_LAI/scripts directory in github

#things you add to the script and can edit when you run, since these will change

#put input of data here as an argument to be used in the script
#vcf_dir="/project/mathilab/gmies/neolithic_selection/1kg_runs/data/full_vcfs/allentoft_snps_vcfs"
#vcf names
#vcf_name="all_samples.chr"
#and within this dir should be vcfs that look like this: $vcf_name${chr}.vcf.gz

# Set input_gen to the value you want, this is used in: ancestry_hmm
input_gen=35 

#set number of admixed individuals, used in: ancestry_hmm
num_admixed=378

#give name of sample file for admixed samples, in data/samples directory outside of results
admixed_samples=asw_samples


#saved num ssa which is $num_admixed * 2
num_admixed_count=$((num_admixed * 2))

#set number of admixed individuals, used in: simplai
num_hg=98
num_neo=48

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




#make directory for ancestry_hmm
mkdir ancestry_hmm
cd ancestry_hmm


#make input files:

#know hg is the minor ancestry here (keeping the files the same)


#make frequency for each ancestry:
module load plink/1.90Beta6.18


#admixed / mneo file:
KEEP_FILE="../input_data/mneo_keep_fam.txt"  # Path to the mneo_keep_fam.txt file

# Initialize output files for the concatenated results
output_neo="neo_combined.frq"
output_hg="hg_combined.frq"
output_bim="combined.bim"

# Clear the output files if they exist
> $output_neo
> $output_hg
> $output_bim

#copy what i need:
cp /project/mathilab/gmies/neolithic_selection/X_chromosome/davy_data/make_data/prunned_cleaned_neo_output.fam full_dataset_chr23_first.fam
cp /project/mathilab/gmies/neolithic_selection/X_chromosome/davy_data/make_data/prunned_cleaned_neo_output.bim full_dataset_chr23.bim
cp /project/mathilab/gmies/neolithic_selection/X_chromosome/davy_data/make_data/prunned_cleaned_neo_output.bed full_dataset_chr23.bed


awk 'BEGIN {FS="\t"; OFS="\t"} {$1="0"; print $0}' full_dataset_chr23_first.fam > full_dataset_chr23.fam


# Process each chromosome
for chr in 23; do
    # Run plink to generate bed file
    #plink --vcf $vcf_dir/$vcf_name${chr}.vcf.gz --const-fid --make-bed --out full_dataset_chr${chr}
    
    # Generate frequency counts for neo dataset
    plink --bfile full_dataset_chr${chr} --keep ../input_data/neo_keep_fam.txt --freq counts --out ${chr}_neo.frq
    
    # Generate frequency counts for hg dataset
    plink --bfile full_dataset_chr${chr} --keep ../input_data/hg_keep_fam.txt --freq counts --out ${chr}_hg.frq
    
    # Concatenate neo frequency files, keeping the header only from the first file
    if [ $chr -eq 1 ]; then
        cat ${chr}_neo.frq.frq.counts > $output_neo
    else
        tail -n +2 ${chr}_neo.frq.frq.counts >> $output_neo
    fi
    
    # Concatenate hg frequency files, keeping the header only from the first file
    if [ $chr -eq 1 ]; then
        cat ${chr}_hg.frq.frq.counts > $output_hg
    else
        tail -n +2 ${chr}_hg.frq.frq.counts >> $output_hg
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
    
    
done



#concat files together:
module load R/4.4
Rscript /project/mathilab/gmies/neolithic_selection/X_chromosome/scripts/ancestry_hmm/concat.R

#combine individuals for each file:

Rscript /project/mathilab/gmies/neolithic_selection/X_chromosome/scripts/ancestry_hmm/combine_indivs.R $num_admixed


#going to delete intermediate files here, so delete the files made from run_chr_mneo.sh
rm *.nosex 
rm *.log
rm *_mneo.frq.counts 


#dont need to add cm to bim, already has it 

#make input file with cm distance and no NAs

Rscript /project/mathilab/gmies/neolithic_selection/X_chromosome/scripts/ancestry_hmm/make_input_file.R ./hg_concat.frq ./neo_concat.frq ./ mneo.counts full_dataset_chr23.bim $num_admixed output_hmm_run_file


#run method:

#run ancestry hmm:
#module load armadillo/12.8.4
#module load openblas
module load Ancestry_HMM/1.0.2


ancestry_hmm -i output_hmm_run_file -s ../input_data/ancestry_hmm_samples.txt -a 2 0.2 0.8 -p 0 0 0.2 -p 1 $input_gen 0.8 -g --ne 10000


#make output file:
mkdir output
cd output

Rscript ../../../scripts/ancestry_hmm/output_average_ancestry_hmm.R ../ output_ancestry_hmm.txt

cp ../output_ancestry_hmm.txt .

cd ..
cd ..


