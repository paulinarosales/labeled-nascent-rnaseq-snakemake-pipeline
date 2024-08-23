log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #
bakr_metaTSV <- snakemake@input[["bakr_metaTSV"]]
feature_ctsTSV <- snakemake@input[["feature_ctsTSV"]]
tx_infoTSV <- snakemake@input[["tx_infoTSV"]]
# genesetTSV <- snakemake@input[["genesetTSV"]]

nascent_ctsTSV <- snakemake@output[["nascent_ctsTSV"]]

min_reads_th <- snakemake@params[["min_reads"]]
min_conv_th <- snakemake@params[["min_conv"]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
conv_tab <- read.table(bakr_metaTSV, stringsAsFactors=FALSE, header = TRUE, sep="\t")
conv_tab <- conv_tab %>% select(c("GF", "T_C"))

feature_cts <- read.table(feature_ctsTSV, stringsAsFactors=FALSE, header = TRUE, sep="\t", comment.char = "#")
feature_cts <- feature_cts %>% select(c(1, 7))
head(feature_cts)
names(feature_cts) <- c("transcript_id", "counts")

tx_info <- read.table(tx_infoTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE)
tx_info <- tx_info %>% select(c("transcript_id", "gene_id"))

# geneset <- read.table(genesetTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE)

# initialize output table
cts_tab <- data.frame(transcript_id = unique(tx_info$transcript_id),
                      total_counts = numeric(length(unique(tx_info$transcript_id))),
                      nascent_counts = numeric(length(unique(tx_info$transcript_id))))
cat("\n")

cat("Assigning total reads from featureCounts...", sep="\n")
cts_tab$total_counts <- feature_cts[match(cts_tab$transcript_id, feature_cts$transcript_id, nomatch = 0), "counts"]
cat("\n")

cat(paste("Filtering out transcript entries with 0 read counts."), sep="\n")
removed <- nrow(cts_tab)
cts_tab <- cts_tab %>% filter(total_counts > 0)
removed <- removed - nrow(cts_tab)
cat(paste("Removed", removed, "entries"), sep="\n")
cat("\n")


# low conversion distribution (for printing only)
tc_counts_dist <- data.frame("T_to_C_conv" = c("0", "1", "2", ">2"),
                             "n_reads" = numeric(4))

tc_counts_dist$n_reads <- c(sum(conv_tab$T_C == 0),
                            sum(conv_tab$T_C == 1),
                            sum(conv_tab$T_C == 2),
                            sum(conv_tab$T_C > 2))

cat(paste0("Filtering out transcript reads with <", min_conv_th, " T-to-C converions nascent reads.\nNascent reads distribution:"), sep="\n")
print(tc_counts_dist)
# count filtering according to min conversions per read threshold
nascent_cts <- conv_tab %>% filter(T_C >= min_conv_th) # remove reads with < min_conv_th T-to-C conversions
cat("\n")

cat("Counting nascent reads...", sep="\n")
nascent_cts <- as.data.frame(table(nascent_cts$GF))
cts_tab$nascent_counts <- nascent_cts[match(cts_tab$transcript_id, nascent_cts$Var1), "Freq"] # ensure right order
cts_tab[is.na(cts_tab)] <- 0  # if transcript_id not found assign 0
cat("\n")

cat(paste("Collapsing transcript counts to gene counts..."), sep="\n") 
cts_tab <- inner_join(tx_info, cts_tab, by="transcript_id")
cts_tab$transcript_id <- NULL

cts_tab <- aggregate(.~gene_id, data = cts_tab, FUN = sum)
cat("\n")

cat("Calculating nascent transcripts fraction...", sep="\n")
cts_tab$nascent_fraction <- round(cts_tab$nascent_counts/cts_tab$total_counts, 2)
cts_tab <- cts_tab[order(cts_tab$nascent_fraction, decreasing = TRUE),]
cat("\n")

cat(paste0("Saving outputs:\n\t- ", nascent_ctsTSV), sep="\n")
write.table(cts_tab, file=nascent_ctsTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")