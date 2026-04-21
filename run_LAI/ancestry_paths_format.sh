#!/bin/bash


#data for ancestralpaths is downloaded from Allentoft et al 2024 VCFs
#then filtered for admixed mneo individuals and snp set 

module load bcftools/1.20 

for chr in {1..22}; do
    echo "Processing chromosome $chr..."

    bcftools query -f '[%AP\t]\n' output_${chr}.neo.vcf.gz | \
    awk '{
        for (i = 1; i <= NF; i++) {
            split($i, a, "|")
            for (j = 1; j <= 2; j++) {
                if (a[j] == "3" || a[j] == "4") {
                    printf "2 "
                } else if (a[j] == "1" || a[j] == "2" || a[j] == "6" || a[j] == "5") {
                    printf "1 "
                } else {
                    printf a[j] " "
                }
            }
        }
        print ""
    }' | gzip > output_ancestrypaths_${chr}.txt.gz

done
