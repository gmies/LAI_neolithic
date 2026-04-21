#!/bin/bash
#1/28/25 to run with sources removed as well! 
#here: /project/mathilab/gmies/neolithic_selection/replication_davy/indep_davy_012825


#need to remake the input files:

#edit from 1kg script to run ancestry hmm

#set number of admixed individuals, used in: ancestry_hmm
num_admixed=367

#make frequency for each ancestry:
module load plink/1.90Beta6.18


#awk '{ if ($3 == "NEO") print $1}' /project/mathilab/gmies/neolithic_selection/qpadm/aadr/extracted_inds/output_eigenstrat_files/neo_output.ind | grep -v -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt > neo_output.ind
#grep -f neo_output.ind /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output.fam | grep -v -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt  | awk '{print $1, $2}' > neo_keep_fam.txt
#awk '{ if ($3 == "HG") print $1}' /project/mathilab/gmies/neolithic_selection/qpadm/aadr/extracted_inds/output_eigenstrat_files/neo_output.ind | grep -v -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt > hg_output.ind
#grep -f hg_output.ind /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output.fam | grep -v -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | awk '{print $1, $2}' > hg_keep_fam.txt


cat \
/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/hg_keep_fam.txt \
/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/neo_keep_fam.txt \
/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/mneo_keep_fam.txt \
| awk '{
    key=$2
    gsub(/\.(AG|SG)$/, "", key)
    gsub(/_(published|enhanced)$/, "", key)
    print key
}' | sort -u > exclude_ids.norm.txt
        
        
        
awk '{ if ($3 == "NEO") print $1}' /project/mathilab/gmies/neolithic_selection/qpadm/aadr/extracted_inds/output_eigenstrat_files/neo_output.ind | grep -v -F -x -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | grep -v -F -f  exclude_ids.norm.txt > neo_output.ind
grep -f neo_output.ind /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output.fam | grep -v -F -f exclude_ids.norm.txt  | grep -v -F -x -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | awk '{print $1, $2}' > neo_keep_fam.txt
awk '{ if ($3 == "HG") print $1}' /project/mathilab/gmies/neolithic_selection/qpadm/aadr/extracted_inds/output_eigenstrat_files/neo_output.ind | grep -v -F -x -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | grep -v -F -f  exclude_ids.norm.txt > hg_output.ind
grep -f hg_output.ind /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output.fam | grep -v -F -f exclude_ids.norm.txt| grep -v -F -x -f  /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | awk '{print $1, $2}' > hg_keep_fam.txt        
        
      
awk '{ if ($3 == "MNEO") print $1}' /project/mathilab/gmies/neolithic_selection/qpadm/aadr/extracted_inds/output_eigenstrat_files/neo_output.ind | grep -v -F -x -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | grep -v -F -f exclude_ids.norm.txt > mneo_output.ind
grep -f mneo_output.ind /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output.fam | grep -v -F -f exclude_ids.norm.txt | grep -v -F -x -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | awk '{print $1, $2}' > mneo_keep_fam.txt  
          
  


plink --bfile /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output --keep neo_keep_fam.txt --extract /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/extract_merged_rsids_davy.txt --freq counts --out neo.frq

plink --bfile /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output --keep hg_keep_fam.txt --extract /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/extract_merged_rsids_davy.txt --freq counts --out hg.frq



#make plink files for mneo indivs

#awk '{ if ($3 == "MNEO") print $1}' /project/mathilab/gmies/neolithic_selection/qpadm/aadr/extracted_inds/output_eigenstrat_files/neo_output.ind | grep -v -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt  > mneo_output.ind
#grep -f mneo_output.ind /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output.fam | grep -v -f /project/mathilab/gmies/neolithic_selection/replication_davy/merge/060225/identical_individuals.txt | awk '{print $1, $2}' > mneo_keep_fam.txt

#ancestry hmm input file:
awk '{print $1, "2"}' mneo_output.ind > ancestry_hmm_samples.txt

#awk '{print $1, "2"}' ../run_ancestry_hmm/ancestry_hmm_samples.txt > ./ancestry_hmm_samples.txt 


#admixed / mneo file:
KEEP_FILE="mneo_keep_fam.txt"  # Path to the mneo_keep_fam.txt file


    #now do for each admixed indiv (mneo)
    # Counter for sequential numbering
    count=1

    # Loop through each line in mneo_keep_fam.txt
    while IFS= read -r line; do
        # Create a temporary keep file for this iteration
        keep_file="${count}_mneo_keep_fam.txt"
        echo "$line" > "$keep_file"

        # Run PLINK command to calculate allele frequencies for the current chromosome
        plink --bfile /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output --keep "$keep_file" --extract /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/extract_merged_rsids_davy.txt --freq counts --out "${count}_mneo"

        # Remove temporary keep file
        rm "$keep_file"

        # Increment the counter
        ((count++))
    done < "$KEEP_FILE"




#concat files together:
module load R/4.4

#combine individuals for each file:

Rscript /project/mathilab/gmies/neolithic_selection/replication_davy/scripts/combine_indivs.R $num_admixed


#going to delete intermediate files here, so delete the files made from run_chr_mneo.sh
rm *.nosex 
rm *.log
rm *_mneo.frq.counts 


#make bim file:
for chr in {1..22}; do
plink --bfile /project/mathilab/gmies/neolithic_selection/admixture/6_12_admixture_wsources/full_admixture_analysis/neo_output --extract /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/extract_merged_rsids_davy.txt --chr $chr --make-bed --out full_dataset_chr${chr}
done

#add cm distance to bim file:

Rscript /project/mathilab/gmies/neolithic_selection/replication_davy/scripts/add_bims_cm.R ./


#make input file with cm distance and no NAs

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/scripts/ancestry_hmm/make_input_file.R ./hg.frq.frq.counts ./neo.frq.frq.counts ./ mneo.counts final_allentoft_bim.bim $num_admixed output_hmm_run_file_dip

       awk 'BEGIN { OFS=" " } {
  for (i = 1; i <= NF; i++) {
    if (i == 1 || i == 2)
      printf "%s", $i;
    else if (i == 7)
      printf "%s", $i;
    else
      printf "%d", $i / 2;

    if (i < NF)
      printf OFS;
    else
      printf "\n";
  }
}' output_hmm_run_file_dip > output_hmm_run_file 

#run method:

#run ancestry hmm:
module load armadillo/12.8.4
module load openblas
module load Ancestry_HMM/1.0.2


ancestry_hmm -i output_hmm_run_file -s ancestry_hmm_samples.txt -a 2 0.19 0.81 -p 0 0 0.19 -p 1 35 0.81 -g --ne 10000

#make output file:
mkdir output
cd output


Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/scripts/ancestry_hmm/output_average_ancestry_hmm.R ../ output_ancestry_hmm.txt

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 ../output_ancestry_hmm.txt


(head -n 1 combined_weighted_z_scores.txt && tail -n +2 combined_weighted_z_scores.txt | sort -k6,6g) | head -n 501 > top_500_sorted_by_pvalue.txt

#plots:

