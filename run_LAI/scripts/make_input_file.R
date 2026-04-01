# This script prepares the input file for ancestry_hmm by combining:
# - Allele counts from two source populations
# - Allele counts for each admixed individual
# - Base pair positions and recombination distances (morgans)
# Inputs:
#   1. Source1 counts file (e.g., hg.frq.counts)
#   2. Source2 counts file (e.g., neo.frq.counts)
#   3. Directory containing admixed sample files
#   4. Sample filename suffix
#   5. Morgan file with recombination distances
#   6. Number of admixed individuals
#   7. Output filename
# Output:
#   - ancestry_hmm-ready input file

library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE)
source1_file <- args[1]
source2_file <- args[2]
input_file_dir <- args[3]
input_file_name <- args[4]
morgan_file_file <- args[5]
num_indivs <- as.numeric(args[6])
output_file <- args[7]

# Read source and Morgan files
source1_counts <- read.table(source1_file, header = TRUE)
source2_counts <- read.table(source2_file, header = TRUE)
morgan_file <- read.table(morgan_file_file, header = FALSE)

# Remove duplicates and filter out SNPs with missing recombination distances
source1_counts <- source1_counts[!duplicated(source1_counts$SNP), ]
source2_counts <- source2_counts[!duplicated(source2_counts$SNP), ]
morgan_file <- morgan_file[!duplicated(morgan_file$V2), ]

morgan_file_saved <- morgan_file
source1_counts <- source1_counts[!is.na(morgan_file$V3), ]
source2_counts <- source2_counts[!is.na(morgan_file$V3), ]
morgan_file <- morgan_file[!is.na(morgan_file$V3), ]

# Merge source populations on CHR, SNP, and alleles
merged_sources <- merge(
  source1_counts, source2_counts,
  by = c("CHR", "SNP", "A1", "A2"),
  all.x = TRUE,
  sort = FALSE
)
merged_sources <- merged_sources[merged_sources$SNP %in% morgan_file$V2, ]

# Add BP and morgans columns
merged_sources$bp <- morgan_file$V4
merged_sources$morgans <- 0

# Compute differences in recombination distances
for (i in 2:nrow(merged_sources)) {
  merged_sources$morgans[i] <- morgan_file$V3[i] - morgan_file$V3[i - 1]
}
merged_sources$morgans[1] <- 0

# Select columns needed for ancestry_hmm
merged_sources_selected <- select(
  merged_sources,
  CHR, bp, C1.x, C2.x, C1.y, C2.y, morgans
)
colnames(merged_sources_selected) <- c("CHR", "bp", "C1.x", "C2.x", "C1.y", "C2.y", "morgans")

# Add allele counts for each admixed individual
for (i in 1:num_indivs) {
  ind_file <- read.table(paste0(input_file_dir, i, input_file_name), header = TRUE)
  ind_file <- ind_file[!duplicated(ind_file$SNP), ]
  ind_file <- ind_file[!is.na(morgan_file_saved$V3), ]
  
  merged_sources_selected[[paste0("sample", i, "_a")]] <- ind_file$C1
  merged_sources_selected[[paste0("sample", i, "_b")]] <- ind_file$C2
}

# Print dimensions
print(dim(merged_sources_selected))

# Write final input file
write.table(merged_sources_selected, output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
