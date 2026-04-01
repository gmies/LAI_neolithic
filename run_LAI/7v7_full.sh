#!/bin/bash


# 7v7_full.sh

# This script runs local and global ancestry inference on the 7v7 sample configuration.
# It executes Ancestry_HMM, RFMix, flare, simpLAI, and MOSAIC pipelines.
# Inputs required:
# - VCF files in ./data/vcf/
# - Sample files (admixed and reference) in ./data/input_files/
# - Map files for flare in ./data/maps/
# Outputs:
# - Results for each method in ./results/{ancestry_hmm, rfmix, flare, simplai, mosaic}/


#-----------------------
# User-defined input
#-----------------------
vcf_dir="./data/vcf"
vcf_name="filt_davy.neo."
input_gen=35 
num_admixed=176
admixed_samples="mneo_samples"

num_admixed_count=$((num_admixed * 2))
num_hg=7
num_neo=7

mkdir -p 7v7
cd 7v7

#-----------------------
# Copy input files for ancestry_hmm, flare
#-----------------------
mkdir -p input_files
cp ./data/input_files/* input_files/

# Reduce reference files to match 7v7 configuration
head -n $num_hg input_files/hg_keep_fam.txt > input_files/hg_keep_fam.txt.tmp && mv input_files/hg_keep_fam.txt.tmp input_files/hg_keep_fam.txt
head -n $((num_hg + num_neo)) input_files/flare.ref.panel > input_files/flare.ref.panel.tmp && mv input_files/flare.ref.panel.tmp input_files/flare.ref.panel

#-----------------------
# Ancestry_HMM
#-----------------------
mkdir -p ancestry_hmm
cd ancestry_hmm

module load plink/1.90Beta6.18

KEEP_FILE="../input_files/mneo_keep_fam.txt"
output_neo="neo_combined.frq"
output_hg="hg_combined.frq"
output_bim="combined.bim"

> $output_neo
> $output_hg
> $output_bim

for chr in {1..22}; do
    plink --vcf $vcf_dir/$vcf_name${chr}.vcf.gz --const-fid --make-bed --out full_dataset_chr${chr}
    plink --bfile full_dataset_chr${chr} --keep ../input_files/neo_keep_fam.txt --freq counts --out ${chr}_neo.frq
    plink --bfile full_dataset_chr${chr} --keep ../input_files/hg_keep_fam.txt --freq counts --out ${chr}_hg.frq

    if [ $chr -eq 1 ]; then
        cat ${chr}_neo.frq.frq > $output_neo
        cat ${chr}_hg.frq.frq > $output_hg
    else
        tail -n +2 ${chr}_neo.frq.frq >> $output_neo
        tail -n +2 ${chr}_hg.frq.frq >> $output_hg
    fi

    cat full_dataset_chr${chr}.bim >> $output_bim

    count=1
    while IFS= read -r line; do
        keep_file="${count}_mneo_keep_fam.txt"
        echo "$line" > "$keep_file"
        plink --bfile "full_dataset_chr${chr}" --keep "$keep_file" --freq counts --out "${chr}_${count}_mneo"
        rm "$keep_file"
        ((count++))
    done < "$KEEP_FILE"
done

module load R/4.2
Rscript ../../scripts/ancestry_hmm/concat.R
Rscript ../../scripts/ancestry_hmm/combine_indivs.R $num_admixed
rm *.nosex *.log *_mneo.frq.counts
Rscript ../../scripts/ancestry_hmm/add_bims_cm.R ./
Rscript ../../scripts/ancestry_hmm/make_input_file.R ./hg_concat.frq ./neo_concat.frq ./ mneo.counts final_allentoft_bim.bim $num_admixed output_hmm_run_file

module load Ancestry_HMM/1.0.2
ancestry_hmm -i output_hmm_run_file -s ../input_files/ancestry_hmm_samples.txt -a 2 0.2 0.8 -p 0 0 0.2 -p 1 $input_gen 0.8 -g --ne 10000

mkdir -p output
cd output
Rscript ../../scripts/ancestry_hmm/output_average_ancestry_hmm.R ../ output_ancestry_hmm.txt
cp ../output_ancestry_hmm.txt .
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output_ancestry_hmm.txt
cd ../..

#-----------------------
# flare
#-----------------------
mkdir -p flare
cd flare
module load bcftools/1.20 htslib gatk

for i in {1..22}; do
    bcftools view -S ../input_files/$admixed_samples -O z -o admixed.chr${i}.vcf.gz $vcf_dir/$vcf_name${i}.vcf.gz
    tabix -p vcf all_samples.chr${i}.vcf.gz
done

for chr in {1..22}; do
    java -jar ./data/flare/flare.jar ref="$vcf_dir/$vcf_name${chr}.vcf.gz" gt="admixed.chr${chr}.vcf.gz" map="./data/maps/plink.chr${chr}.GRCh37.map" ref-panel=../input_files/flare.ref.panel gen=$input_gen min-mac=10 out="${chr}_flare.out"
done

for chr in {1..22}; do
    tabix -p vcf ${chr}_flare.out.anc.vcf.gz
    bcftools query -f '%CHROM\t%POS\t[%AN1]\t[%AN2]\n' ${chr}_flare.out.anc.vcf.gz >> by_indiv_flare_ancestry_output.txt
done

Rscript ../../scripts/flare/make_average.R by_indiv_flare_ancestry_output.txt output_flare.txt

mkdir -p output
mv ../output_flare.txt output/
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output/output_flare.txt
cd ../..

#-----------------------
# mosaic
#-----------------------
cd mosaic
module unload R/4.2
module load R/4.3

Rscript ../../scripts/mosaic/mosaic.R mneo ./ -n $num_admixed -N 10000 --gens $input_gen
Rscript ../../scripts/mosaic/mosaic_output_file.R ./MOSAIC_RESULTS/localanc_mneo_2way_*.RData ./MOSAIC_RESULTS/mneo_2way_*.RData output_mosaic.txt
mkdir -p output
cp output_mosaic.txt output/
cd ..

#-----------------------
# RFMix
#-----------------------
mkdir -p rfmix
cd rfmix

printf "%s " $(yes 0 | head -n 352) > classes.txt
printf "%s " $(yes 1 | head -n 14) >> classes.txt
printf "%s " $(yes 2 | head -n 14) >> classes.txt

for chr in {1..22}; do
    awk '{print substr($0, 1, length($0)-82)}' ./data/rfmix_sep/rfmix_alleles${chr}.txt > rfmix_alleles${chr}.txt
    python2 ./scripts/RunRFMix.py PopPhased rfmix_alleles${chr}.txt classes.txt ./data/rfmix_sep/markerLocationsChr${chr}.txt -G $input_gen -n 5 -o output_rfmix_${chr}
done

Rscript ../../scripts/rfmix/rfmix_output_to_mean.R ./output_rfmix_ ../data/mosaic/rates. output
cp average_ancestry_output.txt output/output_rfmix.txt
mkdir -p output
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output/output_rfmix.txt
cd ../..

#-----------------------
# simpLAI
#-----------------------
mkdir -p simplai
cd simplai
module load R/4.2

num_snps=$(wc -l ../mosaic/mneogenofile.* | grep -v total | awk '{s+=$1} END {print s}')
Rscript ../../scripts/simplai/make_input_simpLAI.R ../mosaic/ $num_snps $num_hg $num_neo $num_admixed input_simplai.gen

num_admixed_count=$((num_admixed * 2))
declare -A genome_lengths=(
    [1]=248956422 [2]=242193529 [3]=198295559 [4]=190214555 [5]=181538259
    [6]=170805979 [7]=159345973 [8]=146364022 [9]=141213431 [10]=135534747
    [11]=135006516 [12]=133851895 [13]=115169878 [14]=107349540 [15]=102531392
    [16]=90354753 [17]=81195210 [18]=78077248 [19]=59128983 [20]=63025520
    [21]=48129895 [22]=50818468
)

for i in {1..22}; do
    mkdir -p ${i}_chr_output
    cd ${i}_chr_output
    head -n 1 ../input_simplai.gen > chr_${i}.gen
    awk -v i="$i" '$1 == i' ../input_simplai.gen >> chr_${i}.gen
    ./data/simpLAI/simpLAI -g chr_${i}.gen --ss1 14 --ss2 14 --ssa $num_admixed_count -l ${genome_lengths[$i]} -s 1e6 -i 5e5 -n 2000 -m 1000 -t 5
    cd ..
done

head -n 1 ../flare/output_flare.txt > full_output_simpLAI.txt
for i in {1..22}; do
    Rscript ../../scripts/simplai/format_simpLAI_output.R ${i}_chr_output/*chr_${i}*fromMin_withSingl.adm $num_admixed_count
    mv simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
    awk -v i="$i" 'NR == 1 {print $0} NR > 1 {$1 = i; print}' ${i}_simpLAI_final_output.txt >> full_output_simpLAI.txt
done

cp full_output_simpLAI.txt output/output_simplai.txt
mkdir -p output
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output/output_simplai.txt
cd ../..
