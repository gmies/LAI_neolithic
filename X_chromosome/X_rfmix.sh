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

#num_hg=48
num_hg=7
num_neo=7



module load plink/1.90Beta6.18

#concat files together:
module load R/4.4


module load bcftools/1.21 
module load htslib

#make input files for mosaic:

#make directory for mosaic
mkdir mosaic_7v7
cd mosaic_7v7

#make input files:


#make pop names:

awk '{print $2, $1}' ../ancestry_hmm_7v7/input_files/flare.ref.panel > sample.names

awk '{print "mneo", $2}' ../ancestry_hmm_7v7/input_files/mneo_keep_fam.txt >> sample.names



#make snp files

awk '{print $1"\t"$2}' ../ancestry_hmm/output/output_ancestry_hmm.txt > keep_davy_snps.txt



# Loop through chromosome numbers 1 to 22

mkdir make_snp_files
cd make_snp_files

chr=23

bcftools view -m2 -M2 -v snps -O v $vcf_dir/$vcf_name${chr}.vcf.gz | bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' > full_no_dist_snpfile.${chr}


#grep -F -f ../keep_davy_snps.txt full_no_dist_snpfile.${chr} > davy_no_dist_snpfile.${chr}

#awk '{print $1"\t"$2"\t""0""\t"$3"\t"$4"\t"$5}' full_no_dist_snpfile.${chr} > ../snpfile.${chr}
awk 'BEGIN{OFS="\t"} {print $2"_"$3"_"$4"_"$5, $2, 0, $3, $4, $5}' full_no_dist_snpfile.${chr} > ../snpfile.${chr}


cd ..



#make recombination map:

Rscript ../scripts/mosaic/make_rates_file.R ./



#make haplotype files:

mkdir haplo_files
cd haplo_files

awk '{print $2}' ../../ancestry_hmm_7v7/input_files/neo_keep_fam.txt > neo_keep.txt
awk '{print $2}' ../../ancestry_hmm_7v7/input_files/hg_keep_fam.txt > hg_keep.txt



# Loop through chromosome numbers 1 to 22

bcftools view -m2 -M2 -v snps ../../flare/admixed.chr${chr}.vcf.gz > filtered_vcf.chr${chr}.vcf

#bcftools query -f '%CHROM\t%POS\t[%GT]\n' filtered_vcf.chr${chr}.vcf > mneo_intermediate_geno.${chr}

#bcftools query -f '%CHROM\t%POS\t[%GT]\n' ../../flare/admixed.chr${chr}.vcf.gz > mneo_intermediate_geno.${chr}


#grep -F -f ../keep_davy_snps.txt mneo_intermediate_geno.${chr} > 2_mneo_intermediate_geno.${chr}

#need to find # columns and get rid of |
#awk '{print $3}' mneo_intermediate_geno.${chr} > 3_mneo_intermediate_geno.${chr}

#sed -i 's/|//g' 3_mneo_intermediate_geno.${chr}



#1/5/26 add double 

bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' filtered_vcf.chr${chr}.vcf \
| awk 'BEGIN{OFS="\t"} {
    for (i=3; i<=NF; i++) {
        if ($i ~ /^[01]$/) {
            $i = $i "|" $i
        }
    }
    print
}' > mneo_intermediate_geno.${chr}


awk '{
    out = ""
    for (i = 3; i <= NF; i++) {
        g = $i
        gsub(/\|/, "", g)
        out = out g
    }
    print out
}' mneo_intermediate_geno.${chr} > 3_mneo_intermediate_geno.${chr}



cp 3_mneo_intermediate_geno.${chr} ../mneogenofile.${chr}


#run with neo and hg:

bcftools view -m2 -M2 -v snps -S neo_keep.txt -O v $vcf_dir/$vcf_name${chr}.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' | awk 'BEGIN{OFS="\t"} {
    for (i=3; i<=NF; i++) {
        if ($i ~ /^[01]$/) {
            $i = $i "|" $i
        }
    }
    print
}'> neo_intermediate_geno.${chr}

#grep -F -f ../keep_davy_snps.txt neo_intermediate_geno.${chr} > 2_neo_intermediate_geno.${chr}

#need to find # columns and get rid of |
#awk '{print $3}' neo_intermediate_geno.${chr} > 3_neo_intermediate_geno.${chr}

#sed -i 's/|//g' 3_neo_intermediate_geno.${chr}

awk '{
    out = ""
    for (i = 3; i <= NF; i++) {
        g = $i
        gsub(/\|/, "", g)
        out = out g
    }
    print out
}' neo_intermediate_geno.${chr} > 3_neo_intermediate_geno.${chr}


cp 3_neo_intermediate_geno.${chr} ../neogenofile.${chr}


#hg:
bcftools view -m2 -M2 -v snps -S hg_keep.txt -O v $vcf_dir/$vcf_name${chr}.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' | awk 'BEGIN{OFS="\t"} {
    for (i=3; i<=NF; i++) {
        if ($i ~ /^[01]$/) {
            $i = $i "|" $i
        }
    }
    print
}' > hg_intermediate_geno.${chr}

#grep -F -f ../keep_davy_snps.txt hg_intermediate_geno.${chr} > 2_hg_intermediate_geno.${chr}

#need to find # columns and get rid of |
#awk '{print $3}' hg_intermediate_geno.${chr} > 3_hg_intermediate_geno.${chr}

#sed -i 's/|//g' 3_hg_intermediate_geno.${chr}

awk '{
    out = ""
    for (i = 3; i <= NF; i++) {
        g = $i
        gsub(/\|/, "", g)
        out = out g
    }
    print out
}' hg_intermediate_geno.${chr} > 3_hg_intermediate_geno.${chr} 


cp 3_hg_intermediate_geno.${chr} ../hggenofile.${chr}


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
printf "%s " $(yes 2 | head -n 14) >> classes.txt


paste -d "" ../mosaic_7v7/mneogenofile.${chr} ../mosaic_7v7/neogenofile.${chr} ../mosaic_7v7/hggenofile.${chr} > rfmix_alleles${chr}.txt

sed -n '3p' ../mosaic_7v7/rates.${chr} | tr ' ' '\n' > markerLocationsChr${chr}.txt

cp -r /project/mathilab/aaw/Mitonuclear2/data/rfmix/PopPhased/ .

sed -i 's/^NA$/0/' markerLocationsChr23.txt


#run method:

module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/3.8



#the number we dont want removed, so if 48 before then -
#90 for 7v3, 82 for 7v7, 94 for 7v1
#awk '{print substr($0, 1, length($0)-82)}' /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/rfmix_alleles${chr}.txt > rfmix_alleles${chr}.txt

python /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/RunRFMix_python2.py PopPhased rfmix_alleles${chr}.txt classes.txt markerLocationsChr${chr}.txt -G $input_gen -n 5 -o output_rfmix_${chr}

#/project/mathilab/aaw/Mitonuclear2/data/rfmix/RunRFMix.py 

#make output file:

Rscript ../scripts/rfmix/rfmix_output_to_mean.R ./output_rfmix_ ../mosaic_7v7/rates. output

cp average_ancestry_output.txt output_rfmix.txt

#cd output
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 output_rfmix.txt


cp output_rfmix_23.0.Viterbi.txt output_rfmix_23.txt.gz
Rscript ../scripts/rfmix/rfmix_global.R ./output_rfmix_


cd ..
cd ..



