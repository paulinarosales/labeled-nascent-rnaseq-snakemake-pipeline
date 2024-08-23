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
convRatesTSV <- snakemake@input[[1]]
outputPDF <- snakemake@output[[1]]

target_mut <- snakemake@params[["base_change"]]
gois <- snakemake@params[["use_genes"]]


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
target_mut <- sub(",", "_", target_mut)
              

conv_tab <- read.table(convRatesTSV, header = TRUE, sep = '\t')
conv_tab <- conv_tab[, c("transcript_name", target_mut)]
conv_tab <- conv_tab %>% filter(transcript_name == gois)
head(conv_tab)
cat("\n")

cat("Plotting...", sep="\n")
rates_p <- ggplot(conv_tab, aes(x=transcript_name, y=T_C, fill=transcript_name)) + 
                geom_bar(position="dodge", stat="identity") +
                geom_text(aes(label = T_C), 
                        angle = 90, position = position_dodge(width = 1), 
                        hjust = -0.2, vjust = 0.5, size = 4) +
                scale_fill_manual(values = pawlette[5:9]) +
                labs(x = "ERCC ID", y = "T-to-C conversion rate [%]") +
                # facet_wrap(~ERCC, nrow = 2) +
                ylim(c(0,0.25)) +
                theme_bw()

cat("\n")



cat(paste0("Saving output:\n\t- ", outputPDF), sep="\n")

pdf(file=outputPDF, width=8, height=6)
                print(rates_p)
        invisible(dev.off())
cat("\n")
cat("DONE!", sep="\n")