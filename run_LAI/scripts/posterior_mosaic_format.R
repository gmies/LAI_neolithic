# Script: posterior_mosaic_format.R
# Purpose: Filter MOSAIC local ancestry calls for high-confidence assignments
#          Converts posterior probabilities to hard calls per SNP per haplotype
#          Outputs compressed per-chromosome files: output_mosaic_<CHR>.txt.gz
# Usage: Rscript mosaic_output_filter.R <RData_path> <g_loc_path> <output_dir>
# Example: Rscript mosaic_output_filter.R ./localanc_file.RData ./g_loc_dir ./output/

library(MOSAIC)
library(data.table)

args <- commandArgs(TRUE)

# Arguments
rdata_file <- args[1]   # Path to MOSAIC RData file containing localanc and g.loc
g_loc_dir   <- args[2]   # Path to directory containing g.loc files
output_dir  <- args[3]   # Directory to save output files

# Load MOSAIC local ancestry RData
load(rdata_file)  # should load `localanc` and `g.loc`
cat("Loaded MOSAIC RData file:", rdata_file, "\n")

# Chromosomes to process
chrnos <- 1:22

# Map local ancestry to SNP positions using MOSAIC function
local_pos <- grid_to_pos(localanc, g_loc_dir, g.loc, chrnos)

# Process each chromosome
for (chr in chrnos) {
  cat("Processing chromosome", chr, "...\n")
  
  # Filter haplotype calls by posterior probability
  # 1 = Neo, 2 = HG
  haplo_calls <- apply(local_pos[[chr]], c(2, 3), function(x) {
    max_prob <- max(x)
    if (max_prob > 0.9) {
      return(which.max(x))  # Keep hard call if probability > 0.9
    } else {
      return(NA)  # Assign NA if max probability <= 0.9
    }
  })
  
  # Transpose to get [sites x haplotypes] as required
  haplo_calls <- t(haplo_calls)
  
  # Create output file path
  out_file <- file.path(output_dir, sprintf("output_mosaic_%d.txt.gz", chr))
  
  # Write haplotype calls to compressed file
  fwrite(haplo_calls, file = out_file, sep = " ", col.names = FALSE,
         row.names = FALSE, quote = FALSE, compress = "gzip")
  
  cat("Saved chromosome", chr, "output to", out_file, "\n")
}

cat("All chromosomes processed successfully.\n")
