#!/bin/bash
# This script runs all Local Ancestry Inference (LAI) methods on a dataset.
# Inputs: VCF files for chromosomes 1-22, sample files, reference panels.
# Outputs: LAI results for ancestry_hmm, RFMix, simpLAI, flare, mosaic, 
#          and sliding-bin Z-scores for each method.
# Usage: Run from the results directory for a specific dataset.

# Configurable inputs
vcf_dir="./data/vcfs"                  # Directory containing input VCFs
vcf_name="filt_davy.neo."              # Prefix for VCF files
input_gen=35                            # Number of generations for ancestry HMM
num_admixed=176                         # Number of admixed individuals
admixed_samples=mneo_samples            # Sample file for admixed individuals

num_admixed_count=$((num_admixed * 2))  # Total alleles
num_hg=48                               # Minor ancestry sample count (HG, CEU)
num_neo=7                               # Major ancestry sample count (NEO, YRI)

# Run ancestry_hmm
cd ancestry_hmm

module load plink/1.90Beta6.18
module load R/4.2

# Prepare input file for ancestry_hmm
Rscript ../../scripts/ancestry_hmm/make_input_file.R \
    ./hg_concat.frq ./neo_concat.frq ./ mneo.counts final_allentoft_bim.bim \
    $num_admixed output_hmm_run_file

# Run ancestry_hmm
module load armadillo/12.8.4
module load openblas
module load Ancestry_HMM/1.0.2

ancestry_hmm -i output_hmm_run_file \
    -s ../input_files/ancestry_hmm_samples.txt \
    -a 2 0.2 0.8 \
    -p 0 0 0.2 \
    -p 1 $input_gen 0.8 -g --ne 10000

# Process ancestry_hmm output
cd output
Rscript ../../../scripts/ancestry_hmm/output_average_ancestry_hmm.R ../ output_ancestry_hmm.txt
cp ../output_ancestry_hmm.txt .

# Compute sliding-bin Z-scores
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 ../output_ancestry_hmm.txt
cd ../..

# Run RFMix
cd rfmix_sep

cp -r ../data/rfmix/PopPhased/ .

module load bcftools/1.20  
module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/3.8
module load R/4.2

for chr in {1..22}; do
    python2 ../scripts/rfmix/RunRFMix.py PopPhased rfmix_alleles${chr}.txt \
        classes.txt markerLocationsChr${chr}.txt -G $input_gen -n 5 \
        -o output_rfmix_${chr}
done

# Process RFMix output
Rscript ../../scripts/rfmix/rfmix_output_to_mean.R ./output_rfmix_ ../mosaic/rates.output
cp average_ancestry_output.txt output/output_rfmix.txt

cd output
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output_rfmix.txt
cd ../..

# Run simpLAI
cd simplai_7_sep

module load R/4.2

num_snps=$(wc -l ../mosaic/mneogenofile.* | grep -v total | awk '{s+=$1} END {print s}')
set -e

declare -A genome_lengths=(
    [1]=248956422 [2]=242193529 [3]=198295559 [4]=190214555
    [5]=181538259 [6]=170805979 [7]=159345973 [8]=146364022
    [9]=141213431 [10]=135534747 [11]=135006516 [12]=133851895
    [13]=115169878 [14]=107349540 [15]=102531392 [16]=90354753
    [17]=81195210 [18]=78077248 [19]=59128983 [20]=63025520
    [21]=48129895 [22]=50818468
)

for i in {1..22}; do
    mkdir -p ${i}_chr_output
    cd ${i}_chr_output

    head -n 1 ../input_simplai.gen > chr_${i}.gen
    awk -v i="$i" '$1 == i' ../input_simplai.gen >> chr_${i}.gen

    ../scripts/simplai/simpLAI -g chr_${i}.gen \
        --ss1 14 --ss2 14 --ssa $num_admixed_count \
        -l ${genome_lengths[$i]} -s 1e6 -i 5e5 -n 2000 -m 1000 -t 5

    cd ..
done

# Combine simpLAI outputs
head -n 1 ../flare/output/output_flare.txt > full_output_simpLAI.txt

for i in {1..22}; do
    Rscript ../../scripts/simplai/format_simpLAI_output.R ${i}_chr_output/*chr_${i}*fromMin_withSingl.adm $num_admixed_count
    mv simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
    awk -v i="$i" 'NR == 1 {print $0} NR > 1 {$1 = i; print}' ${i}_simpLAI_final_output.txt > temp_${i}_simpLAI_final_output.txt
    mv temp_${i}_simpLAI_final_output.txt ${i}_simpLAI_final_output.txt
    awk -v i="$i" '$1 == i' ${i}_simpLAI_final_output.txt >> full_output_simpLAI.txt
done

cp full_output_simpLAI.txt output/output_simplai.txt
cd output
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output_simplai.txt
cd ../..

# Run flare
cd flare

module load bcftools/1.20
module load htslib
module load gatk

for chr in {1..22}; do
    java -jar ../scripts/flare/flare.jar \
        ref="$vcf_dir/$vcf_name${chr}.vcf.gz" \
        gt="admixed.chr${chr}.vcf.gz" \
        map="./data/maps/plink.chr${chr}.GRCh37.map" \
        ref-panel=../input_files/flare.ref.panel \
        gen=$input_gen min-mac=10 out="${chr}_flare.out"
done

for chr in {1..22}; do
    tabix -p vcf ${chr}_flare.out.anc.vcf.gz
    bcftools query -f '%CHROM\t%POS\t[%AN1]\t[%AN2]\n' ${chr}_flare.out.anc.vcf.gz \
        >> by_indiv_flare_ancestry_output.txt
done

Rscript ../../scripts/flare/make_average.R by_indiv_flare_ancestry_output.txt output_flare.txt
cd output
mv ../output_flare.txt .
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output_flare.txt
cd ../..

# Run mosaic
cd mosaic

module load bcftools/1.20
module load R/4.2

Rscript ../../scripts/mosaic/make_rates_file.R ./

module unload R/4.2
module load R/4.3

Rscript ../scripts/mosaic/mosaic.R mneo ./ -n $num_admixed -N 10000 --gens $input_gen
Rscript ../../scripts/mosaic/mosaic_output_file.R ./MOSAIC_RESULTS/localanc_mneo_2way_*.RData \
    ./MOSAIC_RESULTS/mneo_2way_*.RData output_mosaic.txt
cp output_mosaic.txt output/.

cd output
Rscript ../../scripts/sliding_bins_cov_zscores.R 1 output_mosaic.txt
cd ../..
