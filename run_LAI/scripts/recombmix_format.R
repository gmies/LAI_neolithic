#!/usr/bin/env Rscript

library(data.table)

#----------------------------------------------------------
# Input
#----------------------------------------------------------

rfmix_file <- "/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/output/output_rfmix.txt"

recombmix_dir <- "../../recombmix"

rfmix <- fread(rfmix_file)
rfmix[, CHR := as.integer(CHR)]

#----------------------------------------------------------
# Loop over chromosomes
#----------------------------------------------------------

for (chr in 1:22) {

  cat("Processing chromosome", chr, "\n")

  rm_file <- file.path(
    recombmix_dir,
    sprintf("chr%d", chr),
    sprintf("recombmix_chr%d.txt", chr)
  )

  if (!file.exists(rm_file)) {
    cat("Missing:", rm_file, "\n")
    next
  }

  #----------------------------------------------------
  # Parse RecombMix output file
  #----------------------------------------------------

  lines <- readLines(rm_file)
  lines <- lines[nzchar(trimws(lines))]   # drop any blank/whitespace-only lines

  # line 1 = "#Ancestry"
  # following lines = "label\tcode" pairs, until "#Position"
  pos_header_idx <- which(lines == "#Position")

  if (length(pos_header_idx) == 0) {
    warning(sprintf("Chromosome %d: '#Position' header not found, skipping", chr))
    next
  }

  anc_lines <- lines[2:(pos_header_idx - 1)]
  anc_map <- do.call(rbind, strsplit(anc_lines, "\t"))
  colnames(anc_map) <- c("label", "code")
  cat("  Ancestry map:\n")
  print(anc_map)

  # everything after "#Position"
  rest <- lines[(pos_header_idx + 1):length(lines)]

  # positions run until an explicit "#Result" header; haplotype lines follow it
  result_header_idx <- which(rest == "#Result")

  if (length(result_header_idx) == 1) {

    pos_lines <- rest[seq_len(result_header_idx - 1)]
    hap_lines <- rest[(result_header_idx + 1):length(rest)]

  } else {

    # fallback: no "#Result" header found, use the _0/_1 suffix heuristic
    first_field <- sub("\t.*$", "", rest)
    is_hap_line <- grepl("_[01]$", first_field)
    first_hap <- which(is_hap_line)[1]

    if (is.na(first_hap)) {
      warning(sprintf("Chromosome %d: no haplotype lines found, skipping", chr))
      next
    }

    pos_lines <- rest[seq_len(first_hap - 1)]
    hap_lines <- rest[first_hap:length(rest)]
  }

  # parse positions defensively: take only the first tab-delimited field,
  # in case a position line has trailing whitespace/tabs
  positions <- as.integer(sub("\t.*$", "", pos_lines))

  if (anyNA(positions)) {
    bad <- pos_lines[is.na(positions)]
    warning(sprintf(
      "Chromosome %d: %d position line(s) failed to parse as integers, e.g. '%s'",
      chr, length(bad), bad[1]
    ))
  }

  #----------------------------------------------------
  # SNP position table (query points)
  #----------------------------------------------------

  snp_dt <- data.table(POS = positions)
  snp_dt[, start := POS]
  snp_dt[, end := POS]
  setkey(snp_dt, start, end)

  #----------------------------------------------------
  # Per-haplotype interval join
  #----------------------------------------------------

  hap_ids <- character(length(hap_lines))
  anc_list <- vector("list", length(hap_lines))

  for (i in seq_along(hap_lines)) {

    fields <- strsplit(hap_lines[i], "\t")[[1]]
    hap_ids[i] <- fields[1]

    vals <- as.integer(fields[-1])

    if (length(vals) %% 3 != 0) {
      warning(sprintf(
        "Chromosome %d, haplotype %s: segment fields not a multiple of 3, skipping haplotype",
        chr, hap_ids[i]
      ))
      anc_list[[i]] <- rep(NA_integer_, nrow(snp_dt))
      next
    }

    m <- matrix(vals, ncol = 3, byrow = TRUE)
    seg_dt <- data.table(start = m[, 1], end = m[, 2], label = m[, 3])
    setkey(seg_dt, start, end)

    matched <- foverlaps(snp_dt, seg_dt, type = "within", nomatch = NA)

    if (nrow(matched) != nrow(snp_dt)) {
      warning(sprintf(
        "Chromosome %d, haplotype %s: matched %d/%d SNPs",
        chr, hap_ids[i], nrow(matched), nrow(snp_dt)
      ))
    }

    n_na <- sum(is.na(matched$label))
    if (n_na > 0) {
      warning(sprintf(
        "Chromosome %d, haplotype %s: %d SNP(s) fell outside any segment (NA)",
        chr, hap_ids[i], n_na
      ))
    }

    anc_list[[i]] <- matched$label
  }

  anc <- as.data.table(anc_list)
  setnames(anc, hap_ids)

  #----------------------------------------------------
  # Recode: neo (0) -> 1, hg (1) -> 2
  #----------------------------------------------------

  anc <- as.data.table(lapply(anc, function(x) {
    x[x == 0] <- -1L
    x[x == 1] <-  2L
    x[x == -1] <- 1L
    x
  }))

  #----------------------------------------------------
  # Sanity check against rfmix POS list for this chromosome
  #----------------------------------------------------

  rfmix_pos <- rfmix[CHR == chr, POS]

  if (length(positions) != length(rfmix_pos) || !all(positions == rfmix_pos)) {
    warning(sprintf(
      "Chromosome %d: recombmix positions (n=%d) do not match rfmix POS (n=%d)",
      chr, length(positions), length(rfmix_pos)
    ))
  }

  #----------------------------------------------------
  # Write output
  #----------------------------------------------------

  outfile <- sprintf("output_recombmix_%d.txt.gz", chr)

  fwrite(
    anc,
    file = outfile,
    sep = " ",
    col.names = FALSE,
    row.names = FALSE,
    compress = "gzip"
  )

  cat("Wrote", outfile, "\n")
}
