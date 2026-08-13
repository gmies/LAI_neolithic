#7/14/26 to filter gnomix on neolithic data for posterior probabilities
# (analogous to posterior_filtering.R, adapted for gnomix's query_results.fb format)

library(dplyr)
library(data.table)

args <- commandArgs(TRUE)

# arg 1: sample name list (same file used for the RFMix script)
#        ex: /project/mathilab/gmies/.../GA/simplai_7/simplai_global_ancestry.txt
simplai_file <- args[1]

# arg 2: path to the gnomix output base directory, containing chr<N>_output/ subfolders
#        ex: ../../../gnomix   (adjust depending on where you run this from)
gnomix_dir <- args[2]

if (is.na(simplai_file) || is.na(gnomix_dir)) {
  stop("Usage: Rscript gnomix_posterior_filtering.R <simplai_global_ancestry_file> <gnomix_output_dir>")
}

# Read in the sample name list
simplai_data <- read.table(simplai_file, header = TRUE, sep = " ")
sample_names <- simplai_data$sample

all_chromosomes_data <- list()   # per-SNP summary across all haplotypes, one entry per chromosome
all_calls_data <- list()         # per-haplotype calls (SNP rows x haplotype columns), one entry per chromosome

for (chr in 1:22) {

  cat("Processing chromosome", chr, "\n")

  fb_file <- file.path(gnomix_dir, sprintf("chr%d_output", chr), "query_results.fb")

  if (!file.exists(fb_file)) {
    cat("Missing:", fb_file, "\n")
    next
  }

  # line 1 is a comment ("#reference_panel_population: neo hg"), line 2 is the real header
  fb <- fread(fb_file, skip = 1, header = TRUE, sep = "\t")

  # first 4 columns are metadata: chromosome, physical position, genetic_position, genetic_marker_index
  pos_col <- fb[[2]]

  # identify the neo probability column for each haplotype
  # (columns look like "R3:::hap1:::neo", "R3:::hap1:::hg", "R3:::hap2:::neo", ...)
  neo_cols <- grep(":::neo$", colnames(fb), value = TRUE)

  if (length(neo_cols) == 0) {
    warning(sprintf("Chromosome %d: no ':::neo' columns found, skipping", chr))
    next
  }

  neo_probs <- fb[, ..neo_cols]

  # threshold calls: neo prob > 0.9 -> 1 (neo/farmer); neo prob < 0.1 (i.e. hg prob > 0.9) -> 2 (hg); else NA
  calls <- as.data.table(lapply(neo_probs, function(x) {
    out <- rep(NA_integer_, length(x))
    out[x > 0.9] <- 1L
    out[x < 0.1] <- 2L
    out
  }))

  # rename columns to haplotype IDs, e.g. "R3:::hap1:::neo" -> "R3:::hap1"
  setnames(calls, sub(":::neo$", "", neo_cols))

  #------------------------------------------------------
  # per-SNP summary across all haplotypes (like average_ancestry_output.txt)
  #------------------------------------------------------

  chr_summary <- data.table(
    CHR = chr,
    POS = pos_col,
    HGall = rowSums(calls == 2, na.rm = TRUE) /
      (rowSums(calls == 1, na.rm = TRUE) + rowSums(calls == 2, na.rm = TRUE))
  )

  all_chromosomes_data[[chr]] <- chr_summary

  #------------------------------------------------------
  # store haplotype calls to compute by-individual global ancestry later
  #------------------------------------------------------

  all_calls_data[[chr]] <- calls
}

#--------------------------------------------------------
# finish by SNP
#--------------------------------------------------------

final_output <- rbindlist(all_chromosomes_data, use.names = TRUE, fill = TRUE)

write.table(
  final_output,
  file = "average_ancestry_output.txt",
  quote = FALSE,
  row.names = FALSE
)

cat("Wrote average_ancestry_output.txt\n")

#--------------------------------------------------------
# finish by individual
#--------------------------------------------------------

all_data <- rbindlist(all_calls_data, use.names = TRUE, fill = TRUE)

# haplotype column names look like "R3:::hap1", "R3:::hap2"; group by sample ID
hap_cols <- colnames(all_data)
sample_ids <- sub(":::hap[0-9]+$", "", hap_cols)

if (!all(sample_ids %in% sample_names) ) {
  missing_ids <- setdiff(unique(sample_ids), sample_names)
  if (length(missing_ids) > 0) {
    warning(sprintf(
      "The following sample IDs found in gnomix output are not in the sample name list: %s",
      paste(missing_ids, collapse = ", ")
    ))
  }
}

# fraction of hg calls (2) out of all called (non-NA) haplotype calls, per haplotype column
frac_hg <- sapply(hap_cols, function(col) {
  x <- all_data[[col]]
  sum(x == 2, na.rm = TRUE) / sum(!is.na(x))
})

# average the two haplotypes per individual
global_ancestry_by_sample <- tapply(frac_hg, sample_ids, mean, na.rm = TRUE)

# order to match sample_names, and check for completeness
missing_samples <- setdiff(sample_names, names(global_ancestry_by_sample))
if (length(missing_samples) > 0) {
  warning(sprintf(
    "The following samples from the sample list were not found in gnomix output: %s",
    paste(missing_samples, collapse = ", ")
  ))
}

results_df <- data.frame(
  sample = sample_names,
  global_ancestry = global_ancestry_by_sample[sample_names]
)

write.table(
  results_df,
  file = "gnomix_global_ancestry.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Wrote gnomix_global_ancestry.txt\n")
