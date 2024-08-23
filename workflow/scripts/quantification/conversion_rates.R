
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(data.table)
        library(readr)
        library(dplyr)
        })

# ---------- Snakemake parsing ---------- #
collapsedTSV <- snakemake@input[["collapsedTSV"]]
transcriptBED <- snakemake@input[["transcriptBED"]]
globalTSV <- snakemake@output[["globalTSV"]]
txTSV <- snakemake@output[["txTSV"]]

# ---------- Functions ---------- #
computeTxRates <- function(tx_rates_tab, ref_base) {
        cat(paste("Calculating conversion rates for mutations of", ref_base, "as reference base..."), sep="\n")
        base_content <- paste0("n", ref_base)
        regex <- paste0("^", ref_base, "_")
        n_rates <- collapsed_conv %>% select(matches(regex)) # select all combinations of mutations with same original base mutated 
                                                                                # (e.g.  T_C, T_A, T_G for T as original)
        # n_rates$Ncontent <- rowSums(n_rates)
        n_rates <- (n_rates/collapsed_conv[[base_content]]) * 100 # number of mutations of a given base / total number of the original base detected (e.g. T_C/nT)
        print(head(n_rates))
        # n_rates$Ncontent <- NULL
        n_rates <-  round(n_rates, 2)
        return(n_rates)
        cat("\n")
}

# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
collapsed_conv <- read.table(collapsedTSV, header = TRUE, sep = "\t")
# collapsed_conv <- collapsed_conv %>% select(-ends_with("_N"))

tx_info <- read.table(transcriptBED, header = TRUE, sep = "\t")
rownames(tx_info) <- tx_info$transcript_id

tx_rates_tab <- tx_info[collapsed_conv$GF,]
global_rates <- structure(numeric(12), names=c("T_A", "C_A", "G_A",                                                 
                                                "A_T", "C_T", "G_T", 
                                                "A_C", "T_C", "G_C", 
                                                "A_G", "T_G", "C_G"))
cat("\n\n")
cat("Starting with GLOBAL CONVERSION RATES:", sep="\n")
cat("Collapsing mutation counts to global conversion rates...", sep="\n")
global_cts <- collapsed_conv %>% select(c("nA", "nC", "nT", "nG",
                                        "T_A", "C_A", "G_A", "N_A",                                                 
                                        "A_T", "C_T", "G_T", "N_T", 
                                        "A_C", "T_C", "G_C", "N_C", 
                                        "A_G", "T_G", "C_G", "N_G"))

global_cts <- colSums(global_cts)
cat("Global nucleotide content:", sep="\n")
print(global_cts[c("nA", "nC", "nT", "nG")])
cat("\n")

cat("Global mutation counts:", sep="\n")
print(global_cts[c("T_A", "C_A", "G_A", "N_A",                                                 
                   "A_T", "C_T", "G_T", "N_T", 
                   "A_C", "T_C", "G_C", "N_C", 
                   "A_G", "T_G", "C_G", "N_G")])

for(mut in names(global_rates)){
        base_content <- paste0("n", substr(mut, nchar(mut), nchar(mut)))
        global_rates[mut] <- (global_cts[mut]/global_cts[[base_content]]) * 100 # number of mutations of a given base / total number of the original base detected (e.g. T_C/nT)
}
cat("\n")

cat("Global conversion rates:", sep="\n")
global_rates <- data.frame(as.list(round(global_rates, 5)))
print(global_rates)
cat("\n\n")


cat("Starting with INDIVIDUAL RATES:", sep="\n")
cat("Merging transcript counts assignment with conversion counts...", sep="\n")
ref_code <- c("A", "T", "G", "C") # ignore Ns
for(ref_base in ref_code){
        n_rates_tab <- computeTxRates(tx_rates_tab, ref_base)
        tx_rates_tab <- cbind(tx_rates_tab, n_rates_tab)
}
cat("\n")

cat("Removing entries with no conversions...", sep="\n")
print(nrow(tx_rates_tab))
tx_rates_tab <- tx_rates_tab[rowSums(tx_rates_tab[,(ncol(tx_rates_tab)-15):ncol(tx_rates_tab)]) > 0,] # select only cols with the conversion rates
print(nrow(tx_rates_tab))
cat("\n")

cat(paste0("Saving outputs:\n\t- ", globalTSV, "\n\t- ", txTSV), sep="\n")
write.table(global_rates, file = globalTSV, row.names= FALSE, quote = FALSE, sep = "\t")
write.table(tx_rates_tab, file = txTSV, row.names= FALSE, quote = FALSE, sep = "\t")
cat("\n")

cat("DONE!", sep="\n")