#!/bin/bash
#BSUB -J "rfmix_extract[1-22]%6"
#BSUB -n 2
#BSUB -M 8000
#BSUB -W 12:00
#BSUB -o logs/rfmix_chr_%J_%I.out
#BSUB -e logs/rfmix_chr_%J_%I.err

set -euo pipefail

mkdir -p logs

module load bcftools/1.21
module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/3.8
module load R/4.4

chr=${LSB_JOBINDEX}
VCF="../akbari_chr${chr}.vcf.gz"

echo "=== START chr${chr} ==="

if [ ! -f "$VCF" ]; then
    echo "Missing VCF for chr${chr}"
    exit 1
fi

# -------------------------------
# CLEAN OLD FILES (prevents reuse)
# -------------------------------
rm -f ordered_chr${chr}.vcf.gz*
rm -f alleles_chr${chr}.txt
rm -f positions_chr${chr}.txt
rm -f markerLocationsChr${chr}.txt
rm -f rates.${chr}

# -------------------------------
# STEP 1: CREATE CLEAN VCF
# -------------------------------
echo "Creating VCF chr${chr}"

bcftools norm -d any "$VCF" -Oz -o tmp_chr${chr}.vcf.gz
bcftools view -S samples_rfmix_order.txt -Oz -o ordered_chr${chr}.vcf.gz tmp_chr${chr}.vcf.gz
rm tmp_chr${chr}.vcf.gz

bcftools index -f ordered_chr${chr}.vcf.gz

# sanity check
if [ ! -s ordered_chr${chr}.vcf.gz ]; then
    echo "VCF creation failed chr${chr}"
    exit 1
fi

# -------------------------------
# STEP 2: POSITIONS
# -------------------------------
echo "Extracting positions chr${chr}"

bcftools query -f '%POS\n' ordered_chr${chr}.vcf.gz > positions_chr${chr}.txt

# -------------------------------
# STEP 3: ALLELES
# -------------------------------
echo "Extracting alleles chr${chr}"

bcftools query -f '[%GT\t]\n' ordered_chr${chr}.vcf.gz | \
sed 's/[|\/]/ /g' | \
awk '{
    for(i=1;i<=NF;i++){
        printf "%s", substr($i,1,1)
        printf "%s", substr($i,2,1)
    }
    printf "\n"
}' > alleles_chr${chr}.txt

# -------------------------------
# STEP 4: CHECK SNP COUNTS
# -------------------------------
n1=$(wc -l < positions_chr${chr}.txt)
n2=$(wc -l < alleles_chr${chr}.txt)

echo "SNPs chr${chr}: positions=$n1 alleles=$n2"

if [ "$n1" -ne "$n2" ]; then
    echo "ERROR: SNP mismatch chr${chr}"
    exit 1
fi

# -------------------------------
# STEP 5: RATES
# -------------------------------
echo "Generating rates chr${chr}"

Rscript make_rates_file_from_vcf.R ${chr} positions_chr${chr}.txt ./

if [ ! -s rates.${chr} ]; then
    echo "Rates failed chr${chr}"
    exit 1
fi

# -------------------------------
# STEP 6: MARKER LOCATIONS
# -------------------------------
sed -n '3p' rates.${chr} | tr ' ' '\n' > markerLocationsChr${chr}.txt

n3=$(wc -l < markerLocationsChr${chr}.txt)

echo "Marker count chr${chr}: $n3"

if [ "$n3" -ne "$n1" ]; then
    echo "ERROR: markerLocations mismatch chr${chr}"
    exit 1
fi

# -------------------------------
# STEP 7: RFMIX
# -------------------------------
echo "Running RFMix chr${chr}"

python /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/RunRFMix_python2.py \
    PopPhased \
    alleles_chr${chr}.txt \
    classes.txt \
    markerLocationsChr${chr}.txt \
    -G 35 \
    -o output_rfmix_${chr}

echo "=== DONE chr${chr} ==="
