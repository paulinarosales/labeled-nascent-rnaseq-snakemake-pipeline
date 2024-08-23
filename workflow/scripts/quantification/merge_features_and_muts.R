# Adapted from: https://github.com/simonlabcode/bam2bakR/blob/main/workflow/scripts/merge_features_and_muts.R

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(data.table)
        library(readr)
        library(dplyr)
        library(tidyverse)
        })

# ---------- Snakemake parsing ---------- #
featureCountsTSV <- snakemake@input[["featureCountsTSV"]]
conversionCountsTSV <- snakemake@input[["conversionCountsTSV"]]
metaTSV <- snakemake@output[["metaTSV"]]
collapsedTSV <- snakemake@output[["collapsedTSV"]]


cat("Reading input data...", sep="\n")
feature_cts <- fread(featureCountsTSV)
colnames(feature_cts) <- c("qname", "status", "nhits", "GF")
conversion_cts <- fread(conversionCountsTSV, check.names = FALSE)
cat("\n")

cat("Merging transcript counts assignment with conversion counts...", sep="\n")
cat(paste0("\tUsing files:\n\t\t- ", featureCountsTSV, "\n\t\t- ", conversionCountsTSV), sep="\n")
head(conversion_cts)
cat("Removing unaligned reads...", sep="\n")
feature_cts <- feature_cts[ nhits > 0 , c("qname", "GF")] # remove no alignments
feature_cts[, GF := gsub(",", "+", GF)]
head(feature_cts)
conversion_cts <- feature_cts[conversion_cts, on = .(qname)]
conversion_cts <- conversion_cts %>% drop_na(GF)
# colnames: qname	GF	nA	nC	nT	nG	rname	FR	sj	TA	CA	GA	V12	AT	CT	GT	NT	AC	TC	GC	NC	AG	TG	CG	NG	AN	TN	CN	GN	NN	gmutloc	tp

# conversion_cts$qname <- NULL
# conversion_cts$rname <- NULL
# conversion_cts$sj <- NULL
# conversion_cts$gmutloc <- NULL
# conversion_cts$tp <- NULL
cat("\n")

collapsed_cts <- conversion_cts %>% select(-c(qname, FR, rname, sj, gmutloc, tp))
head(collapsed_cts)
collapsed_cts <- aggregate(.~GF, data = collapsed_cts, FUN = sum)
head(collapsed_cts)
cat(paste0("Saving outputs:\n\t- ", metaTSV), sep="\n")
write.table(conversion_cts, file = metaTSV, row.names= FALSE, quote = FALSE, sep = "\t")
write.table(collapsed_cts, file = collapsedTSV, row.names= FALSE, quote = FALSE, sep = "\t")
cat("\n")

cat("DONE!", sep="\n")