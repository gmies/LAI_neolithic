#!/bin/bash

#7/7/26

#run recombmix 


mkdir shared_LAI_inputs
cd shared_LAI_inputs


# add in:
vcf_dir="/project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/data"
#vcf names
vcf_name="chrX.phased."

module load bcftools/1.21

module load htslib/1.21 

module load boost/1.74.0  

######################################################
# sample lists
######################################################

awk '{print $2}' ../ancestry_hmm_7v7/input_files/neo_keep_fam.txt > reference_samples.txt
awk '{print $2}' ../ancestry_hmm_7v7/input_files/hg_keep_fam.txt >> reference_samples.txt

######################################################
# GNomix sample map
######################################################

echo -e "#Sample\tPanel" > gnomix.smap

awk '{print $2"\tneo"}' ../ancestry_hmm_7v7/input_files/neo_keep_fam.txt >> gnomix.smap
awk '{print $2"\thg"}' ../ancestry_hmm_7v7/input_files/hg_keep_fam.txt >> gnomix.smap

######################################################
# RecombMix labels
######################################################

echo -e "#Sample_ID\tPopulation_Label" > recombmix.labels

awk '{print $2"\tneo"}' ../ancestry_hmm_7v7/input_files/neo_keep_fam.txt >> recombmix.labels
awk '{print $2"\thg"}' ../ancestry_hmm_7v7/input_files/hg_keep_fam.txt >> recombmix.labels

######################################################
# chromosome files
######################################################

chr=23


echo "Preparing chromosome ${chr}"

bcftools view \
    -S reference_samples.txt \
    -Oz \
    -o reference.chr${chr}.vcf.gz \
    $vcf_dir/$vcf_name${chr}.vcf.gz

tabix -p vcf reference.chr${chr}.vcf.gz

awk 'BEGIN{OFS="\t";print "#chromosome","position","Genetic_Map(cM)"} \
{print $1,$4,$3}' \
/project/mathilab/gmies/neolithic_selection/flare/maps/plink.chrX.GRCh37.map \
> gnomix_map.chr${chr}.tsv


cd ..

