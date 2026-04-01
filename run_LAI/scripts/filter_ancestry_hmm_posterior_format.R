# Script to process RFMix/ancestry_hmm posterior files per individual and chromosome
# Outputs a file per chromosome with ancestry states coded as "2 2", "1 2", "1 1", or "NA NA"
# Started 4/1/26

library(tidyverse)

# Paths to input/output
individuals_file <- "../data/mneo_samples"             # Text file listing individual IDs
posterior_dir <- "../data/ancestry_hmm"                # Directory containing .posterior.gz files
output_dir <- "../data/output_ancestryhmm"             # Directory to write per-chromosome output

# Read the list of individuals
individuals <- readLines(individuals_file)

# Loop through each chromosome 1-22
for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  # Initialize a list to store ancestry data for all individuals for this chromosome
  chr_data <- list()
  
  # Loop through each individual
  for (ind in individuals) {
    
    # Construct path to this individual's posterior file
    posterior_file <- file.path(posterior_dir, paste0(ind, ".posterior.gz"))
    
    # Check if file exists
    if (file.exists(posterior_file)) {
      
      # Read posterior probabilities
      posterior_data <- read.table(posterior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      
      # Keep only rows for the current chromosome
      chr_rows <- posterior_data[posterior_data$chrom == chr, ]
      
      # Extract the last three columns (ancestry probabilities)
      posterior_values <- chr_rows[, 3:5]
      
      # Determine ancestry states per SNP based on max probability
      ancestry_states <- apply(posterior_values, 1, function(row) {
        max_val <- max(row)
        
        if (max_val <= 0.9) {
          return("NA NA")       # Low confidence
        } else if (max_val == row[1]) {
          return("2 2")         # Ancestry state 2 homozygous
        } else if (max_val == row[2]) {
          return("1 2")         # Heterozygous
        } else {
          return("1 1")         # Ancestry state 1 homozygous
        }
      })
      
      # Store this individual's ancestry states
      chr_data[[ind]] <- ancestry_states
      
    } else {
      cat("Warning: File not found:", posterior_file, "\n")
    }
  }
  
  # Combine all individuals' ancestry states for this chromosome into a matrix
  chr_data_combined <- do.call(cbind, chr_data)
  
  # Write combined ancestry states to a gzipped file
  output_file <- file.path(output_dir, paste0("output_ancestryhmm_", chr, ".txt.gz"))
  write.table(chr_data_combined, output_file, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat("Finished processing chromosome", chr, "\n")
}
