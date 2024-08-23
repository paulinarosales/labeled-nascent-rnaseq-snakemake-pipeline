log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
	library(gridExtra)
        library(ggplot2)
        library(ggvenn)
        library(GGally)
})

# ---------- Snakemake parsing ---------- #
cts_sitesTSV <- snakemake@input[["cts_sitesTSV"]]
ref_geneListTSV <- snakemake@input[["ref_geneListTSV"]]
ref_ctsTSV <- snakemake@input[["ref_ctsTSV"]]
genesetTSV <- snakemake@input[["genesetTSV"]]

vennPDF <- snakemake@output[["vennPDF"]]
pairsPDF <- snakemake@output[["pairsPDF"]]

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

sites_t <- read.table(cts_sitesTSV, header = TRUE, sep = "\t", row.names = )
ref_genes <- read.table(ref_geneListTSV, header = TRUE, sep = "\t")
ref_cts <- read.table(ref_ctsTSV, header = TRUE, sep = "\t")
geneset <- read.table(genesetTSV, header = TRUE, sep = "\t")

ref_cts <- as.data.frame(table(ref_cts$gene))
colnames(ref_cts) <- c("gene", "mut_count")
head(ref_cts)
cat("\n")

# ----- Venn diagrams
cat("Identifying overlapping genes with m6A sites...", sep="\n")

overlap_genes <- list(Ref = ref_genes$gene, 
                      Sample = sites_t$gene)

cat("Ploting Venn Diagram...", sep="\n")
p_venn <- ggvenn(overlap_genes, fill_color = pawlette[1:2],
                    stroke_size = 0.5, set_name_size = 4) +
                    labs(title = basename(dirname(cts_sitesTSV))) +
                    theme(plot.title = element_text(hjust = 0.5, face="bold")) 

pdf(file=vennPDF, width=8, height=6)
                print(p_venn)
        invisible(dev.off())
cat(paste0("Plot saved at:\n\t- ", vennPDF), sep="\n")
cat("\n")

# ----- Correlation scatter plots

cat("Merging with sites counts from reference...", sep="\n")
gene_names <- unique(geneset$gene_name)
merged_cts <- data.frame("Ref"=numeric(length(gene_names)),
                         "Sample"=numeric(length(gene_names)),
                         row.names = gene_names)
head(merged_cts)
merged_cts$Ref <- ref_cts[match(rownames(merged_cts), ref_cts$gene), "mut_count"]

merged_cts$Sample <- sites_t[match(rownames(merged_cts), sites_t$gene), "mut_count"]
head(merged_cts)
cat("Ploting correlations...", sep="\n")

# p_corr <- ggpairs(merged_cts, #columnLabels = gsub('_', ' ', colnames(vsd_rep1), fixed = T), 
#                   labeller = label_wrap_gen(10), 
#                   title = basename(dirname(cts_sitesTSV)),
#                   diag=list(continuous='blank'), 
#                   lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.5)))

p_corr <- ggpairs(merged_cts, title = basename(dirname(cts_sitesTSV)))


pdf(file=pairsPDF, width=8, height=6)
                print(p_corr)
        invisible(dev.off())
cat(paste0("Plot saved at:\n\t- ", pairsPDF), sep="\n")
cat("\n")



cat("\n")
cat("DONE!", sep="\n")



	# curPlot = ggplot(plotTab, aes(x=class,y=values,fill=highlight,col=highlight)) + 
        # stat_boxplot(geom ='errorbar') + 
        # geom_boxplot(outlier.shape = NA,lwd=0.8,fatten=2) + 
        # facet_grid(~group, scales="free", space="free") + 
        # xlab("") + ylab("Mutation rate per UTR base [%]") +
	# scale_fill_manual(values=c("white","white")) + 
        # scale_color_manual(values=c("black", "red")) + 
        # theme(axis.ticks.x = element_blank(), legend.position = "none") + 
        # coord_cartesian(ylim=c(0, ymax))

	# plotList[[length(plotList)+1]] <- curPlot + ggtitle(rates$sample[i])


