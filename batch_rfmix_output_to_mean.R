#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(parallel)

args <- commandArgs(TRUE)

# arg1 = prefix to RFMix files (e.g. "./output_rfmix_")
dir_file_path <- args[1]

# number of cores (set manually or detect)
n_cores <- detectCores() - 1

process_chr <- function(i) {

  cat("Processing chr", i, "\n")

  # --- RFMix file ---
  rfmix_file <- paste0(dir_file_path, i, ".0.Viterbi.txt")
  chr_txt <- fread(rfmix_file, header = FALSE)

  # --- VCF POS extraction (fast) ---
  vcf_file <- paste0("ordered_chr", i, ".vcf.gz")

  chr_map <- fread(
    cmd = paste("bcftools query -f '%POS\\n'", vcf_file),
    col.names = "POS"
  )

  # --- sanity check ---
  if (nrow(chr_txt) != nrow(chr_map)) {
    stop(paste("Mismatch in chr", i,
               "| RFMix rows:", nrow(chr_txt),
               "| VCF rows:", nrow(chr_map)))
  }

  # --- ancestry calculation ---
  num_2 <- rowSums(chr_txt == 2, na.rm = TRUE)
  num_1 <- rowSums(chr_txt == 1, na.rm = TRUE)

  data.frame(
    CHR = i,
    POS = chr_map$POS,
    HGall = num_2 / (num_1 + num_2)
  )
}

# --- run in parallel ---
results_list <- mclapply(1:22, process_chr, mc.cores = n_cores)

# --- combine in correct order ---
final_output <- bind_rows(results_list)


#remove duplicate snps:
# --- deduplicate by CHR + POS ---
final_output <- final_output %>%
  group_by(CHR, POS) %>%
  summarise(HGall = mean(HGall, na.rm=TRUE), .groups='drop') %>%
  arrange(CHR, POS)




# --- write once ---
write.table(final_output,
            file = "average_ancestry_output.txt",
            quote = FALSE,
            row.names = FALSE)


cat("Done! Output written to average_ancestry_output.txt\n")
