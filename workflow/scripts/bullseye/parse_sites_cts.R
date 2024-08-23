log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
        library(tidyr)
})

# ---------- Snakemake parsing ---------- #

racTSV <- snakemake@input[[1]]

parsedTSV <- snakemake@output[["parsedTSV"]]
geneListTSV <- snakemake@output[["geneListTSV"]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
sites_tab <- read.table(racTSV, header = TRUE,  sep = '\t', comment.char = "")
colnames(sites_tab)[1] <- "chr"
cat("\n")
head(sites_tab)  
cat("Parsing input table...", sep="\n")
parse_tab <- separate(sites_tab, 
                    col=gene_region_etc, into=c("gene", "region", "mut", "mut_count", 
                                                "control_edit_ratio_out", "type"), sep='\\|')

head(parse_tab)      

parse_tab <- parse_tab %>% select(-c("chr", "start", "end", "strand"))
parse_tab$mut_count <- sub("mut=", "", parse_tab$mut_count)
head(parse_tab)  

gene_tab <- parse_tab %>% select(c("gene", "mut_count"))
head(gene_tab)
gene_tab$mut_count <- as.numeric(gene_tab$mut_count)
gene_tab <-  aggregate(.~gene, data = gene_tab, FUN = sum)
head(gene_tab)
# colnames(gene_tab) <- c("gene", "mut_count")
cat("\n")

cat(paste0("Saving outputs:\n\t- ", parsedTSV), sep="\n")
write.table(parse_tab, file=parsedTSV, sep="\t", quote=FALSE, row.names=FALSE)
write.table(gene_tab, file=geneListTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")