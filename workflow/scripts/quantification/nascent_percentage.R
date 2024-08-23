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

nascent_percentTSV <- snakemake@output[[1]]

# ---------- Main code ---------- #

sample_names <- basename(dirname(cts_files))
nascent_content <- data.frame(sample = sample_names,
                              lib_percent = numeric(length(sample_names)),
                              genes_percent = numeric(length(sample_names)),
                              row.names = sample_names)


cat(paste("Calculating nascent reads content for all samples..."), sep="\n")
# Merge rates from all samples
for(i in 1:length(cts_files)) {
    sample_file <- cts_files[i]
    sample_name <- sample_names[i]

    cat(paste0(sample_name, ":"), sep="\n")
    cat(paste("\tReading counts table..."), sep="\n")
    sample_cts <- read.table(sample_file, header = TRUE, sep = '\t')

    cat(paste("\tCalculating total library percentage corresponding to nascent reads..."), sep="\n")
    nascent_lib <- (sum(sample_cts$nascent_counts) / sum(sample_cts$total_counts)) * 100
    nascent_content[sample_name, "lib_percent"] <- round(nascent_lib, 1)
    cat("\n")

    cat(paste("\tFiltering out gene entries with any reads"), sep="\n")
    keep <- rowSums(sample_cts[,c("total_counts", "nascent_counts")]) > 0 # remove entries with any reads
    sample_cts <- sample_cts[keep,]
    cat("\n")
    
    cat(paste("\tIdentifying nascent genes..."), sep="\n")
    positive_genes <- sum(sample_cts$nascent_counts > 0, na.rm = TRUE)
    nascent_genes <- (positive_genes / nrow(sample_cts)) * 100
    nascent_content[sample_name, "genes_percent"] <- round(nascent_genes, 1)
    cat(paste("\t\tFound", positive_genes, " genes with at least 1 nascent read."), sep="\n")
    print(nrow(sample_cts))
    cat("\n")
}
cat("\n")

print(nascent_content)

cat(paste0("Saving output:\n\t- ", nascent_percentTSV), sep="\n")
write.table(nascent_content, file=nascent_percentTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")