# concat.R
#-------------------------------------------------------
# Purpose: Concatenate per-chromosome allele frequency counts
#          separately for hg and neo populations
# Input: Files named {chr}_hg.frq.frq.counts and {chr}_neo.frq.frq.counts (chr 1–22)
# Output: Two combined files: hg_concat.frq and neo_concat.frq
#-------------------------------------------------------

hg_concat <- NULL
neo_concat <- NULL

# Loop over chromosomes 1-22
for (chr in 1:22) {
  
  # Construct file names
  hg_file <- paste0(chr, "_hg.frq.frq.counts")
  neo_file <- paste0(chr, "_neo.frq.frq.counts")
  
  # Read in the chromosome data
  hg_data <- read.table(hg_file, header = TRUE, stringsAsFactors = FALSE)
  neo_data <- read.table(neo_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Concatenate chromosome data
  hg_concat <- if (is.null(hg_concat)) hg_data else rbind(hg_concat, hg_data)
  neo_concat <- if (is.null(neo_concat)) neo_data else rbind(neo_concat, neo_data)
}

# Write combined output files
write.table(hg_concat, "hg_concat.frq", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(neo_concat, "neo_concat.frq", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
