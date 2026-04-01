#!/usr/bin/env Rscript

# filter_ancestry_hmm_format.R
# ----------------------
# Reads per-individual ancestry HMM posterior probability files (.posterior.gz)
# Converts probabilities into discrete ancestry calls for each SNP per chromosome
# Outputs per-chromosome haplotype call files (tab-separated, no header, compressed)
#
# Usage:
# Rscript ./scripts/filter_ancestry_hmm_format.R <individuals_file> <posterior_dir> <output_dir>
# Example:
# Rscript ./scripts/filter_ancestry_hmm_format.R ./data/mneo_samples ./data/ancestry_hmm ./output

library(data.table)

args <- commandArgs(TRUE)

if(length(args) != 3) {
  stop("Usage: Rscript ancestryhmm_to_haplo.R <individuals_file> <posterior_dir> <output_dir>")
}

individuals_file <- args[1]   # Text file with sample names, one per line
posterior_dir    <- args[2]   # Directory containing .posterior.gz files
output_dir       <- args[3]   # Directory to save output files

# Read the list of individuals
individuals <- readLines(individuals_file)

# Process chromosomes 1–22
for(chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  chr_data <- list()  # Will store haplotype calls for all individuals
  
  for(ind in individuals) {
    posterior_file <- file.path(posterior_dir, paste0(ind, ".posterior.gz"))
    
    if(!file.exists(posterior_file)) {
      cat("Warning: File not found:", posterior_file, "\n")
      next
    }
    
    # Read posterior probabilities
    posterior_data <- fread(posterior_file, header = TRUE, sep = "\t", data.table = FALSE)
    
    # Filter rows for the current chromosome
    chr_rows <- posterior_data[posterior_data$chrom == chr, ]
    
    # Extract the last three columns (ancestry probabilities)
    posterior_values <- chr_rows[, (ncol(chr_rows)-2):ncol(chr_rows)]
    
    # Convert to discrete ancestry calls
    ancestry_states <- apply(posterior_values, 1, function(row) {
      max_val <- max(row)
      if(max_val <= 0.9) {
        return("NA NA")        # Low confidence, mark as missing
      } else if(max_val == row[1]) {
        return("2 2")          # HG/HG
      } else if(max_val == row[2]) {
        return("1 2")          # Neo/HG
      } else {
        return("1 1")          # Neo/Neo
      }
    })
    
    chr_data[[ind]] <- ancestry_states
  }
  
  # Combine all individuals for this chromosome
  chr_data_combined <- do.call(cbind, chr_data)
  
  # Write output as tab-separated, no header, compressed
  output_file <- file.path(output_dir, paste0("output_ancestryhmm_", chr, ".txt.gz"))
  fwrite(chr_data_combined, file = output_file, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE, compress = "gzip")
  
  cat("Finished processing chromosome", chr, "\n")
}
