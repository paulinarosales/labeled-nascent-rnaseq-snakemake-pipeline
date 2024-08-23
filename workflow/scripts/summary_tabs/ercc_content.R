log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #

idxstatsTSV <- snakemake@input[[1]]
idx_percentTSV <- snakemake@output[["idx_percentTSV"]]
ercc_summaryTSV <- snakemake@output[["ercc_summaryTSV"]]

cat("Reading input data...", sep="\n")
stats_t <- read.table(idxstatsTSV, header = FALSE, sep = '\t')
colnames(stats_t) <- c("seq_name", "length", "n_mapped", "n_unmapped")
cat("\n")

cat(paste("Removing sequences with 0 mapped reads-segments..."), sep="\n")
stats_t <- stats_t %>% filter(n_mapped > 0)
cat("\n")

cat(paste("Calculating mapped reads percentage per sequence..."), sep="\n")
stats_t$mapped_percent <- (stats_t$n_mapped/sum(stats_t$n_mapped))*100
print(sum(stats_t$mapped_percent))
cat("\n")

cat(paste("Calculating ERCC reads content..."), sep="\n")
ercc_t <- stats_t[startsWith(stats_t$seq_name, "ERCC-"),]
head(ercc_t)

summary_t <- data.frame("Total_mapped_reads" = sum(ercc_t$n_mapped), 
                        "Percentage_mapped_reads" = sum(ercc_t$mapped_percent))

head(summary_t)
cat("\n")


cat(paste0("Saving outputs:\n\t- ", idx_percentTSV, "\n\t- ", ercc_summaryTSV), sep="\n")
write.table(stats_t, file = idx_percentTSV, row.names= FALSE, quote = FALSE, sep = "\t")
write.table(summary_t, file = ercc_summaryTSV, row.names= FALSE, quote = FALSE, sep = "\t")
cat("\n")

cat("DONE!", sep="\n")