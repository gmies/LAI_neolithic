#!/bin/bash
set -euo pipefail

module load plink/2.0-20210505
module load plink/1.90Beta6.18
module load admixture

#    PATHS   
VCF_DIR="/project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/data"
VCF_NAME="chrX.phased.23.vcf.gz"



#make input files:
#mkdir input_files
cp /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/input_files/* input_files/.

#num of hg we want 
head -n 3 input_files/hg_keep_fam.txt > input_files/hg_keep_fam.tmp && mv input_files/hg_keep_fam.tmp input_files/hg_keep_fam.txt


#sum of hg and farmers we want 
head -n 10 input_files/flare.ref.panel > input_files/flare.ref.panel.tmp && mv input_files/flare.ref.panel.tmp input_files/flare.ref.panel

cat input_files/mneo_keep_fam.txt input_files/neo_keep_fam.txt input_files/hg_keep_fam.txt > keep_targets_plus_sources.txt

cat input_files/neo_keep_fam.txt input_files/hg_keep_fam.txt > 1kg_ceu.yri.fid_iid.txt

KEEP="keep_targets_plus_sources.txt"   # targets + refs
OUTDIR="./output_chrX_admixture"

mkdir -p ${OUTDIR}
cd ${OUTDIR}


mv ../*.txt .

#    STEP 1   
echo "=== Make chrX PLINK dataset ==="

plink \
  --vcf ${VCF_DIR}/${VCF_NAME} \
  --keep ${KEEP} \
  --const-fid \
  --make-bed \
  --out chrX.all

# NOTE:
# males will be auto-encoded as diploid (0/0 or 1/1)

#    STEP 2   
echo "=== LD pruning (chrX) ==="

plink \
  --bfile chrX.all \
  --indep-pairwise 200 50 0.1 \
  --out chrX.pruned

plink \
  --bfile chrX.all \
  --extract chrX.pruned.prune.in \
  --make-bed \
  --out chrX.pruned

#    STEP 3   
echo "=== Run ADMIXTURE ==="

K=2
admixture --cv chrX.pruned.bed ${K} | tee admixture_chrX_k${K}.out




# map IID → row number in .fam
awk '
NR==FNR {row[$2]=FNR; next}
($2 in row) {print row[$2]}
' chrX.pruned.fam 1kg_ceu.yri.fid_iid.txt > source_rows.txt

# remove sources from Q
awk '
NR==FNR {bad[$1]=1; next}
!(FNR in bad)
' source_rows.txt chrX.pruned.${K}.Q > targets_only.${K}.Q




awk '{print $2}' chrX.pruned.fam > sample_ids.txt

paste sample_ids.txt chrX.pruned.${K}.Q > chrX_admixture_with_ids.txt

