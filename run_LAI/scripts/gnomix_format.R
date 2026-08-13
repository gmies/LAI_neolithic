#!/usr/bin/env Rscript

library(data.table)

#----------------------------------------------------------
# Input
#----------------------------------------------------------


rfmix_file <- "/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/output/output_rfmix.txt"

gnomix_dir <- "../../gnomix"

rfmix <- fread(rfmix_file)
rfmix[, CHR := as.integer(CHR)]

#----------------------------------------------------------
# Loop over chromosomes
#----------------------------------------------------------

for(chr in 1:22){

    cat("Processing chromosome", chr, "\n")

    msp_file <- file.path(
        gnomix_dir,
        sprintf("chr%d_output/query_results.msp", chr)
    )

    if(!file.exists(msp_file)){
        cat("Missing:", msp_file, "\n")
        next
    }

    #------------------------------------------------------
    # Read MSP
    #------------------------------------------------------

    msp <- fread(
        msp_file,
        skip = 1,          # skip "#Subpopulation..."
        header = TRUE
    )

    # Clean annoying column names
    names(msp) <- gsub(" ", "_", names(msp))
    names(msp)[1:6] <- c("chr","spos","epos","sgpos","egpos","nsnps")

    # Interval table
    msp[, start := spos]
    msp[, end   := epos]
    setkey(msp, start, end)

    #------------------------------------------------------
    # SNP positions
    #------------------------------------------------------

    snps <- data.table(
        POS = rfmix[CHR == chr, POS]
    )

    if(nrow(snps)==0){
        cat("No SNPs\n")
        next
    }

    snps[, start := POS]
    snps[, end := POS]
    setkey(snps, start, end)

    #------------------------------------------------------
    # Interval join
    #------------------------------------------------------

    matched <- foverlaps(
        snps,
        msp,
        type="within",
        nomatch=0L
    )

    if(nrow(matched)!=nrow(snps)){
        warning(sprintf(
            "Chromosome %d: matched %d/%d SNPs",
            chr,
            nrow(matched),
            nrow(snps)
        ))
    }

    #------------------------------------------------------
    # Extract ancestry calls
    #------------------------------------------------------


print(names(matched))   # run once to confirm exact names in your data

meta_cols <- c("POS", "i.start", "i.end", "start", "end",
               "chr", "spos", "epos", "sgpos", "egpos", "nsnps")

anc_cols <- setdiff(names(matched), meta_cols)
anc <- matched[, ..anc_cols]

#    anc <- matched[, -(1:10), with=FALSE]

    # Convert coding:
    # neo:0 -> 2
    # hg :1 -> 1

anc <- as.data.table(lapply(anc, function(x){
    x[x == 0] <- -1L   # temporary value
    x[x == 1] <-  2L
    x[x == -1] <- 1L
    x
}))

    #------------------------------------------------------
    # Write output
    #------------------------------------------------------

    outfile <- sprintf("output_gnomix_%d.txt.gz", chr)

    fwrite(
        anc,
        file=outfile,
        sep=" ",
        col.names=FALSE,
        row.names=FALSE,
        compress="gzip"
    )

    cat("Wrote", outfile, "\n")
}
