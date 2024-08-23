log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
        library(tidyr)
	library(gridExtra)
        library(ggplot2)
})

# ---------- Snakemake parsing ---------- #
convFiles <- snakemake@input[["convFiles"]]
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]
all_globalRatesTSV <- snakemake@output[["all_globalRatesTSV"]]
barplotPDF <- snakemake@output[["barplotPDF"]]

# ---------- Color palette ---------- #
pawlette <- c("#e18727", # orange
		"#0072b5", # dark blue
		"#20854e", # dark green
		"#bc3c29", # red
		"#7876b1", # purple
		"#cd6090", # pink
		"#ffdb24", # yellow
		"#6ca6cd", # light blue
		"#cd6090", # pink
		"#698b22", # light green
		"#8b5a2b") # brown

# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
sample_man <- read.table(sample_manifestTSV, header = TRUE, sep = "\t")
# rownames(sample_man) <- basename(dirname(convFiles)) #paste(sample_man$Sample_type, sample_man$Treatment, "Bio-rep", sample_man$Bio_rep, sep="_")
sample_man$files <- convFiles
cat("\n")

cat("Merging conversion rates tables for all samples...", sep="\n")
rates_t <- data.frame("T_A" = numeric(nrow(sample_man)), "C_A" = numeric(nrow(sample_man)), "G_A" = numeric(nrow(sample_man)),                                                 
                      "A_T" = numeric(nrow(sample_man)), "C_T" = numeric(nrow(sample_man)), "G_T" = numeric(nrow(sample_man)), 
                      "A_C" = numeric(nrow(sample_man)), "T_C" = numeric(nrow(sample_man)), "G_C" = numeric(nrow(sample_man)), 
                      "A_G" = numeric(nrow(sample_man)), "T_G" = numeric(nrow(sample_man)), "C_G" = numeric(nrow(sample_man)),
                      "Sample_type" = sample_man$Sample_type, 
                      "Treatment" = sample_man$Treatment,
                      "Sample" = basename(dirname(convFiles)))



for(i in 1:nrow(sample_man)){
        conv_rates <- read.table(sample_man$files[i], header = TRUE, sep = "\t")
        rates_t[i, colnames(conv_rates)] <- conv_rates
}


pivot_rates_t <- rates_t %>% pivot_longer(cols = 1:12, names_to = "mut", values_to = "rate")
head(pivot_rates_t)
pivot_rates_t$rate <- round(pivot_rates_t$rate, 2)
head(pivot_rates_t)
# pivot_rates_t$mut <- sub("_", "-to-", pivot_rates_t$mut) # chage mutation label from N_N --> N-to-N for clarity


rates_t <- rates_t %>% select(c("Sample", 1:12)) # select only rate to save plot
cat("\n")

cat("Plotting...", sep="\n")
rates_p <- ggplot(pivot_rates_t, aes(x=mut, y=rate, fill=Sample_type, alpha=Treatment)) + 
                geom_bar(position="dodge", stat="identity") +
                geom_text(aes(label = rate, group=Treatment), 
                        angle = 90, position = position_dodge(width = 1), 
                        hjust = -0.2, vjust = 0.5, size = 2, alpha = 1) +
                scale_fill_manual(values = pawlette[1:length(unique(pivot_rates_t$Sample_type))]) +
                scale_alpha_manual(values = c(1, 0.4)) +
                labs(x = "Mutation", y = "Global conversion rate [%]") +
                facet_wrap(~Sample_type, ncol=1) +
                ylim(c(0,1)) +
                theme_bw()

cat("\n")
cat(paste0("Saving outputs:\n\t- ", all_globalRatesTSV, "\n\t-", barplotPDF), sep="\n")

write.table(rates_t, file = all_globalRatesTSV, row.names= FALSE, quote = FALSE, sep = "\t")
pdf(file=barplotPDF, width=8, height=6)
                print(rates_p)
        invisible(dev.off())
cat("\n")
cat("DONE!", sep="\n")