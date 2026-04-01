#!/bin/bash

#11/12/25 update to run on differing sample sizes of sources and compare GA correlations (and LA as well)

#script to run all lai methods on a dataset as an input 
#start this script in the results section of a directory


#script 6/11/25 to rerun with gens 35 for all methods 
module load bcftools/1.20
module load htslib


#things you add to the script and can edit when you run, since these will change

#put input of data here as an argument to be used in the script
vcf_dir="/project/mathilab/gmies/neolithic_selection/allentoft_data/ancestry_prop/data/031825_filtering/davy_snps"
#vcf names
vcf_name="filt_davy.neo."
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

#num_hg = 48
num_hg=1
num_neo=7

mkdir 7v1
cd 7v1

mkdir input_files
cp /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/* input_files/.

#num of hg we want 
head -n 1 input_files/hg_keep_fam.txt > input_files/hg_keep_fam.tmp && mv input_files/hg_keep_fam.tmp input_files/hg_keep_fam.txt


#sum of hg and farmers we want 
head -n 8 input_files/flare.ref.panel > input_files/flare.ref.panel.tmp && mv input_files/flare.ref.panel.tmp input_files/flare.ref.panel



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

# Process each chromosome
for chr in {1..22}; do
    # Run plink to generate bed file
    plink --vcf $vcf_dir/$vcf_name${chr}.vcf.gz --const-fid --make-bed --out full_dataset_chr${chr}
    
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
    
    
done

#concat files together:
module load R/4.2

Rscript ../../scripts/ancestry_hmm/concat.R


#combine individuals for each file:

Rscript ../../scripts/ancestry_hmm/combine_indivs.R $num_admixed

#going to delete intermediate files here, so delete the files made from run_chr_mneo.sh
rm *.nosex 
rm *.log
rm *_mneo.frq.counts 



#add cm distance to bim file:

Rscript ../../scripts/ancestry_hmm/add_bims_cm.R ./




#know hg is the minor ancestry here (keeping the files the same)


#make frequency for each ancestry:
module load plink/1.90Beta6.18
module load R/4.2

#make input file with cm distance and no NAs

Rscript ../../scripts/ancestry_hmm/make_input_file.R ./hg_concat.frq ./neo_concat.frq ./ mneo.counts final_allentoft_bim.bim $num_admixed output_hmm_run_file




#run method:

#run ancestry hmm:
module load armadillo/12.8.4
module load openblas
module load Ancestry_HMM/1.0.2


ancestry_hmm -i output_hmm_run_file -s ../input_files/ancestry_hmm_samples.txt -a 2 0.2 0.8 -p 0 0 0.2 -p 1 $input_gen 0.8 -g --ne 10000




#make output file:
mkdir output
cd output

Rscript ../../../scripts/ancestry_hmm/output_average_ancestry_hmm.R ../ output_ancestry_hmm.txt

cp ../output_ancestry_hmm.txt .


#run sliding bins z scores:
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 ../output_ancestry_hmm.txt

cd ..
cd ..



#make input files:
module load R/4.2


#flare:

mkdir flare
cd flare

#make input files:

# Loop through chromosomes 1-22
for i in {1..22}
do

bcftools view -S /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/data/samples/$admixed_samples -O z -o admixed.chr${i}.vcf.gz $vcf_dir/$vcf_name${i}.vcf.gz

tabix -p vcf all_samples.chr${i}.vcf.gz 
    
done


#make just ASW vcfs

module load bcftools/1.20
module load htslib


#run method:

#run flare by chromosome: 

module load gatk
# Loop through chromosome numbers 1 to 22
for chr in {1..22}
do
java -jar /project/mathilab/gmies/neolithic_selection/flare/flare.jar ref="$vcf_dir/$vcf_name${chr}.vcf.gz" gt="admixed.chr${chr}.vcf.gz" map="/project/mathilab/gmies/neolithic_selection/flare/maps/plink.chr${chr}.GRCh37.map" ref-panel=../input_files/flare.ref.panel gen=$input_gen  min-mac=10 out="${chr}_flare.out"

done

#had removed these previously: excludemarkers=/project/mathilab/gmies/neolithic_selection/flare/remove_snps/remove.variants.list


#make output file:

#combine results and make LAI by pos:


module load bcftools/1.20 
module load htslib
# Loop through chromosome numbers 1 to 22
for chr in {1..22}
do

tabix -p vcf ${chr}_flare.out.anc.vcf.gz

bcftools query -f '%CHROM\t%POS\t[%AN1]\t[%AN2]\n' ${chr}_flare.out.anc.vcf.gz >> by_indiv_flare_ancestry_output.txt

done

#then make output in r:

module load R/4.2
Rscript ../../scripts/flare/make_average.R by_indiv_flare_ancestry_output.txt output_flare.txt


cd output

mv ../output_flare.txt .

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 output_flare.txt

cd ..



cd ..



#make input files:
module load R/4.2





#make input files for mosaic:

#make directory for mosaic
mkdir mosaic
cd mosaic

#make input files:


#make pop names:

awk '{print $2, $1}' ../input_files/flare.ref.panel > sample.names

awk '{print "mneo", $2}' ../input_files/mneo_keep_fam.txt >> sample.names



#make snp files

awk '{print $1"\t"$2}' ../ancestry_hmm/output/output_ancestry_hmm.txt > keep_davy_snps.txt


module load bcftools/1.20 

# Loop through chromosome numbers 1 to 22

mkdir make_snp_files
cd make_snp_files

for chr in {1..22}

do

bcftools view -m2 -M2 -v snps -O v $vcf_dir/$vcf_name${chr}.vcf.gz | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' > full_no_dist_snpfile.${chr}


grep -F -f ../keep_davy_snps.txt full_no_dist_snpfile.${chr} > davy_no_dist_snpfile.${chr}

awk '{print $1"\t"$2"\t""0""\t"$3"\t"$4"\t"$5}' davy_no_dist_snpfile.${chr} > ../snpfile.${chr}

done

cd ..



#make recombination map:

module load R/4.2

Rscript ../../scripts/mosaic/make_rates_file.R ./



#make haplotype files:

mkdir haplo_files
cd haplo_files

awk '{print $2}' ../../input_files/neo_keep_fam.txt > neo_keep.txt
awk '{print $2}' ../../input_files/hg_keep_fam.txt > hg_keep.txt


module load bcftools/1.20 

# Loop through chromosome numbers 1 to 22

for chr in {1..22}

do

bcftools view -m2 -M2 -v snps ../../flare/admixed.chr${chr}.vcf.gz > filtered_vcf.chr${chr}.vcf

bcftools query -f '%CHROM\t%POS\t[%GT]\n' filtered_vcf.chr${chr}.vcf > mneo_intermediate_geno.${chr}

#bcftools query -f '%CHROM\t%POS\t[%GT]\n' ../../flare/admixed.chr${chr}.vcf.gz > mneo_intermediate_geno.${chr}


grep -F -f ../keep_davy_snps.txt mneo_intermediate_geno.${chr} > 2_mneo_intermediate_geno.${chr}

#need to find # columns and get rid of |
awk '{print $3}' 2_mneo_intermediate_geno.${chr} > 3_mneo_intermediate_geno.${chr}

sed -i 's/|//g' 3_mneo_intermediate_geno.${chr}

cp 3_mneo_intermediate_geno.${chr} ../mneogenofile.${chr}


#run with neo and hg:

bcftools view -m2 -M2 -v snps -S neo_keep.txt -O v $vcf_dir/$vcf_name${chr}.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' > neo_intermediate_geno.${chr}

grep -F -f ../keep_davy_snps.txt neo_intermediate_geno.${chr} > 2_neo_intermediate_geno.${chr}

#need to find # columns and get rid of |
awk '{print $3}' 2_neo_intermediate_geno.${chr} > 3_neo_intermediate_geno.${chr}

sed -i 's/|//g' 3_neo_intermediate_geno.${chr}

cp 3_neo_intermediate_geno.${chr} ../neogenofile.${chr}


#hg:
bcftools view -m2 -M2 -v snps -S hg_keep.txt -O v $vcf_dir/$vcf_name${chr}.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' > hg_intermediate_geno.${chr}

grep -F -f ../keep_davy_snps.txt hg_intermediate_geno.${chr} > 2_hg_intermediate_geno.${chr}

#need to find # columns and get rid of |
awk '{print $3}' 2_hg_intermediate_geno.${chr} > 3_hg_intermediate_geno.${chr}

sed -i 's/|//g' 3_hg_intermediate_geno.${chr}

cp 3_hg_intermediate_geno.${chr} ../hggenofile.${chr}


done

cd ..
cd ..



#make directory for rfmix
mkdir rfmix
cd rfmix

#make input files:

printf "%s " $(yes 0 | head -n 352) > classes.txt 

#num farmers *2
printf "%s " $(yes 1 | head -n 14) >> classes.txt 

#num hg *2
printf "%s " $(yes 2 | head -n 2) >> classes.txt




cp -r /project/mathilab/aaw/Mitonuclear2/data/rfmix/PopPhased/ .



#run method:

module load bcftools/1.20  
module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/3.8
module load R/4.2


for chr in {1..22}
do

#the number we want 
awk '{print substr($0, 1, length($0)-94)}' /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/rfmix_alleles${chr}.txt > rfmix_alleles${chr}.txt

python2 /project/mathilab/aaw/Mitonuclear2/data/rfmix/RunRFMix.py PopPhased rfmix_alleles${chr}.txt classes.txt /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/markerLocationsChr${chr}.txt -G $input_gen -n 5 -o output_rfmix_${chr}

done


#make output file:

Rscript ../../scripts/rfmix/rfmix_output_to_mean.R ./output_rfmix_ ../mosaic/rates. output

cp average_ancestry_output.txt output/output_rfmix.txt

cd output
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 output_rfmix.txt

cd ..
cd ..




#and run simplai 7v7



#and run simplai 7v7
#i think just edit to grab the first 7 of the hg (ceu)
#make directory for simplai
mkdir simplai
cd simplai

#make input files:
module load R/4.2

#this is edited by number of individuals, would need to change for more individuals ***
#and need num snps in mosaic output ***
num_snps=$(wc -l ../mosaic/mneogenofile.* | grep -v total | awk '{s+=$1} END {print s}')

Rscript ../../scripts/simplai/make_input_simpLAI.R ../mosaic/ $num_snps $num_hg $num_neo $num_admixed input_simplai.gen

#run method:

#saved num ssa which is $num_admixed * 2
num_admixed_count=$((num_admixed * 2))

# Ensure script stops on error
set -e

# Define the genome lengths for each chromosome (you can modify these based on actual genome sizes)
declare -A genome_lengths=(
    [1]=248956422
    [2]=242193529
    [3]=198295559
    [4]=190214555
    [5]=181538259
    [6]=170805979
    [7]=159345973
    [8]=146364022
    [9]=141213431
    [10]=135534747
    [11]=135006516
    [12]=133851895
    [13]=115169878
    [14]=107349540
    [15]=102531392
    [16]=90354753
    [17]=81195210
    [18]=78077248
    [19]=59128983
    [20]=63025520
    [21]=48129895
    [22]=50818468
)

for i in {1..22}; do
    # Create output directory for each chromosome
    mkdir -p ${i}_chr_output
    cd ${i}_chr_output

   # Extract and write the header (first line) to the output file
    head -n 1 ../input_simplai.gen > chr_${i}.gen
    
    # Extract data for the specific chromosome
    awk -v i="$i" '$1 == i' ../input_simplai.gen >> chr_${i}.gen
    
    # Run simpLAI with the extracted data, using the chromosome-specific genome length
    /project/mathilab/gmies/neolithic_selection/simpLAI/simpLAI -g chr_${i}.gen --ss1 2 --ss2 14 --ssa $num_admixed_count -l ${genome_lengths[$i]} -s 1e6 -i 5e5 -n 2000 -m 1000 -t 5
    
    # Go back to the parent directory
    cd ..
done



#make output file:


module load R/4.2


head -n 1 ../flare/output_flare.txt > full_output_simpLAI.txt

for i in {1..22}; do
Rscript ../../scripts/simplai/format_simpLAI_output.R ${i}_chr_output/*chr_${i}*fromMin_withSingl.adm $num_admixed_count

mv simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
awk -v i="$i" 'NR == 1 {print $0} NR > 1 {$1 = i; print}' ${i}_simpLAI_final_output.txt > temp_${i}_simpLAI_final_output.txt
mv temp_${i}_simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
awk -v i="$i" '$1 == i' ${i}_simpLAI_final_output.txt >> full_output_simpLAI.txt

done

cp full_output_simpLAI.txt output/output_simplai.txt

mkdir output
cd output
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 output_simplai.txt

cd ..


cd ..








#run mosaic:
#run method:

cd mosaic

module unload R/4.2
module load R/4.3

Rscript /project/mathilab/gmies/MOSAIC-master/mosaic.R mneo ./ -n $num_admixed -N 10000 --gens $input_gen 


#make output file:


Rscript ../../scripts/mosaic/mosaic_output_file.R ./MOSAIC_RESULTS/localanc_mneo_2way_*.RData ./MOSAIC_RESULTS/mneo_2way_*.RData output_mosaic.txt  

module unload R/4.3

mkdir output
cp output_mosaic.txt output/.


cd ..

























