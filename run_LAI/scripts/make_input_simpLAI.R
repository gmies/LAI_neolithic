# make_input_simpLAI.R
#-------------------------------------------------------
# Purpose: Construct input file for running simpLAI
# Inputs: 
#   1) input_dir: directory containing SNP and genotype files
#   2) num_snps: total number of SNPs
#   3) num_hg: number of hunter-gatherer individuals
#   4) num_neo: number of Neolithic individuals
#   5) num_mneo: number of admixed Neolithic individuals
#   6) output_file_name: desired name for output file
# Output: Tab-separated file containing SNPs and genotype columns for all groups
#-------------------------------------------------------

library(dplyr)

args <- commandArgs(TRUE)
input_dir <- args[1]
num_snps <- as.numeric(args[2])
num_hg <- as.numeric(args[3])
num_neo <- as.numeric(args[4])
num_mneo <- as.numeric(args[5])
output_file_name <- args[6]

#-------------------------------------------------------
# Generate column names for each group (hg, neo, mneo)
#-------------------------------------------------------
generate_column_names <- function(group_prefix, num_individuals) {
  colnames <- c()
  for (i in 1:num_individuals) {
    colnames <- c(colnames, 
                  paste0(group_prefix, "i", i, "_1"), 
                  paste0(group_prefix, "i", i, "_2"))
  }
  return(colnames)
}

hg_columns <- generate_column_names("s1_", num_hg)
neo_columns <- generate_column_names("s2_", num_neo)
mneo_columns <- generate_column_names("sa_", num_mneo)

all_columns <- c("Chrom", "Pos", "Anc_all", "Der_all", hg_columns, neo_columns, mneo_columns)

#-------------------------------------------------------
# Initialize empty data frame for final LAI file
#-------------------------------------------------------
LAI_file <- data.frame(matrix(NA, nrow = num_snps, ncol = length(all_columns)))
colnames(LAI_file) <- all_columns

#-------------------------------------------------------
# Fill in data chromosome by chromosome
#-------------------------------------------------------
saved_snp_start <- 1

for (i in 1:22) {
  # Read SNP file
  snp_file <- read.table(paste0(input_dir, "snpfile.", i), header = FALSE)
  small_lai_file <- data.frame(snp_file$V2, snp_file$V4, snp_file$V5, snp_file$V6)
  
  LAI_file[saved_snp_start:(saved_snp_start + nrow(snp_file) - 1), 1:4] <- small_lai_file
  
  # Read genotypes for each group as fixed-width files
  hg_geno <- read.fwf(paste0(input_dir, "hggenofile.", i), widths = rep(1, num_hg*2), header = FALSE)
  neo_geno <- read.fwf(paste0(input_dir, "neogenofile.", i), widths = rep(1, num_neo*2), header = FALSE)
  mmneo_geno <- read.fwf(paste0(input_dir, "mneogenofile.", i), widths = rep(1, num_mneo*2), header = FALSE)
  
  # Fill in columns for each group
  LAI_file[saved_snp_start:(saved_snp_start + nrow(snp_file) - 1), 5:(4 + num_hg*2)] <- hg_geno
  LAI_file[saved_snp_start:(saved_snp_start + nrow(snp_file) - 1), (5 + num_hg*2):(4 + num_hg*2 + num_neo*2)] <- neo_geno
  LAI_file[saved_snp_start:(saved_snp_start + nrow(snp_file) - 1), (5 + num_hg*2 + num_neo*2):(4 + num_hg*2 + num_neo*2 + num_mneo*2)] <- mmneo_geno
  
  saved_snp_start <- saved_snp_start + nrow(snp_file)
}

#-------------------------------------------------------
# Write output file
#-------------------------------------------------------
write.table(LAI_file, output_file_name, quote = FALSE, row.names = FALSE, col.names = TRUE)
