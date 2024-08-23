log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #

cts_files <- snakemake@input[["cts_files"]]
genesetTSV <- snakemake@input[["genesetTSV"]]

merged_ctsTSV <- snakemake@output[[1]]

merge_col <- snakemake@wildcards[["counts"]]
min_reads_th <- snakemake@params[["min_reads_th"]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
merged_cts <- read.table(genesetTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE) # initialize merged_cts with gene info

sample_names <- basename(dirname(cts_files))
merged_cts[sample_names] <- numeric(nrow(merged_cts))
cat("\n")

cat(paste0("Merging ", merge_col, " for all samples..."), sep="\n")
# Merge rates from all samples
for(i in 1:length(cts_files)) {
  #i = 1
    sample_file <- cts_files[i]
    sample_name <- sample_names[i]
    sample_cts <- read.table(sample_file, header = TRUE, sep = '\t')

    merged_cts[[sample_name]] <- sample_cts[match(merged_cts$gene_id, sample_cts$gene_id), merge_col]
}
cat("\n")

merged_cts[is.na(merged_cts)] <- 0

if(merge_col != "total_counts"){
    min_reads_th <- 0
    keep <- rowSums(merged_cts[,sample_names]) > min_reads_th
}else{
    keep <- rowSums(merged_cts[,sample_names]) >= min_reads_th
}

cat(paste("Filtering out gene entries with <", min_reads_th, "total reads across samples."), sep="\n")
removed <- nrow(merged_cts)
merged_cts <- merged_cts[keep,]
removed <- removed - nrow(merged_cts)
cat(paste("Removed", removed, "gene entries"), sep="\n")
cat("\n")

cat(paste0("Saving outputs:\n\t- ", merged_ctsTSV), sep="\n")
write.table(merged_cts, file=merged_ctsTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")