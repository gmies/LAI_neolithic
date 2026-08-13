#!/bin/bash

#gnomix to rfmix format


module load R/4.4

#!/bin/bash

#gnomix to rfmix format


module load R/4.4

Rscript gnomix_format.R

Rscript ../mosaic_7v7/rfmix_output_to_mean.R ./output_gnomix_ ../mosaic_7v7/rates.

Rscript /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/remove_duplicate_cols.R output_gnomix_23.txt.gz

Rscript ../mosaic_7v7/rfmix_output_to_mean.R ./output_gnomix_23_final ../mosaic_7v7/rates.
Rscript ../mosaic_7v7/rfmix_global.R ./output_gnomix_

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 average_ancestry_output.txt

python  /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "output_gnomix_\${chr}_final.txt.gz" \
  --output-prefix gnomix_chr23 \
  --chr-start 23 \
  --chr-end 23
