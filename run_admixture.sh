#!/bin/bash
set -euo pipefail

module load plink/1.90Beta6.18
mkdir admixture
cd admixture

# ---- Input files ----
MNEO="../input_files/mneo_keep_fam.txt"
NEO="../input_files/neo_keep_fam.txt"
HG="../input_files/hg_keep_fam.txt"

# ---- Output keep file ----
KEEP="keep_fam.txt"

# ---- PLINK/ADMIXTURE input files ----
DIR="/project/mathilab/gmies/neolithic_selection/admixture/allentoft_admixture/imputed/2_20_25"
BASE="${DIR}/prunned_cleaned_neo_output"

FAM="${BASE}.fam"
BED="${BASE}.bed"
BIM="${BASE}.bim"

# ---- Output filtered prefix ----
OUT="prunned_cleaned_neo_output.keep"

# ---- Binaries ----
PLINK="/project/mathilab/bin/plink"
ADMIXTURE="/project/mathilab/bin/admixture"


echo "=== Step 1: Create keep_fam.txt ==="
cat "$MNEO" "$NEO" "$HG" > "$KEEP"
echo "Wrote: $KEEP"


echo "=== Step 2: Running PLINK --keep to make filtered dataset ==="
plink \
  --bfile "$BASE" \
  --keep "$KEEP" \
  --make-bed \
  --out "$OUT"


echo "Filtered files generated:"
ls -lh ${OUT}.bed ${OUT}.bim ${OUT}.fam


# ---- Verify ADMIXTURE input exists ----
if [[ ! -f "${OUT}.bed" ]]; then
    echo "ERROR: ${OUT}.bed was not created. Cannot run ADMIXTURE."
    exit 1
fi


echo "=== Step 3: Running ADMIXTURE (K=2) ==="
echo "Command:"
echo "$ADMIXTURE ${OUT}.bed 2"

$ADMIXTURE "${OUT}.bed" 2

echo "ADMIXTURE complete!"
echo "Generated:"
echo "  ${OUT}.2.Q"
echo "  ${OUT}.2.P"
