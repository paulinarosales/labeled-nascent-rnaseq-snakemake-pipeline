log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
	library(gridExtra)
        library(ggplot2)
})

# ---------- Snakemake parsing ---------- #
conRatesTSV <- snakemake@input[[1]]
target_mut <- snakemake@params[["base_change"]]
outputPDF <- snakemake@output[[1]]

# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
rates_t <- read.table(conRatesTSV, header = TRUE, sep = "\t")

#       chromosome	start	end	
#       strand	transcript_name	
#       gene_id	gene_name       gene_type	
#       A_T	A_C	A_G	
#       T_A	T_C	T_G	
#       G_A	G_T	G_C	
#       C_A	C_T	C_G


target_mut <- sub(",", "_", target_mut)
cat("\n")

rates_t <- rates_t %>% select(c("A_T","A_C","A_G","T_A","T_C","T_G","G_A","G_T","G_C","C_A","C_T","C_G"))
head(rates_t)
quantiles <- lapply(rates_t, function(x) {
				return(quantile(x, na.rm=TRUE, p=0.75) + 1.5 * IQR(x, na.rm=TRUE))
		})

print(quantiles)

ymax <- ceiling(max(unlist(quantiles)))
# ymax <- 3


rates_t <- rbind(# A
		data.frame(class = "A_T", values = rates_t$A_T),
		data.frame(class = "A_C", values = rates_t$A_C),
		data.frame(class = "A_G", values = rates_t$A_G),
                # T
		data.frame(class = "T_A", values = rates_t$T_A),
		data.frame(class = "T_C", values = rates_t$T_C),
		data.frame(class = "T_G", values = rates_t$T_G),
                # G
		data.frame(class = "G_A", values = rates_t$G_A),
		data.frame(class = "G_T", values = rates_t$G_T),
		data.frame(class = "G_C", values = rates_t$G_C),
                # C
		data.frame(class = "C_A", values = rates_t$C_A),
		data.frame(class = "C_T", values = rates_t$C_T),
		data.frame(class = "C_G", values = rates_t$C_G)
)
head(rates_t)

cat(paste0("Identifying target mutation: ", sub("_", ">", target_mut)), sep="\n")
rates_t$highlight <- "no"
rates_t$highlight[rates_t$class == target_mut] <- "yes"
rates_t$class <- sub("_", ">", rates_t$class)

# Divide in original base blocks
rates_t$group <- "A"
rates_t$group[rates_t$class %in% c("T>A","T>C","T>G")] <- "T"
rates_t$group[rates_t$class %in% c("G>A","G>T","G>C")] <- "G"
rates_t$group[rates_t$class %in% c("C>A","C>T","C>G")] <- "C"
head(rates_t)
rates_t <- rates_t[!is.na(rates_t$values),]
cat("\n")

cat("Ploting...", sep="\n")
rates_p <- ggplot(rates_t, aes(x=class,y=values,fill=highlight,col=highlight)) + 
                stat_boxplot(geom ='errorbar') + 
                geom_boxplot(outlier.shape = NA,lwd=0.8,fatten=2) + 
                facet_grid(~group, scales="free", space="free") + 
                xlab("") + 
                ylab("Mutation rate per transcript base [%]") +
		scale_fill_manual(values=c("white","white")) + 
                scale_color_manual(values=c("black", "red")) + 
                theme(axis.ticks.x = element_blank(), legend.position = "none") + 
                coord_cartesian(ylim=c(0, ymax))
cat("\n")
cat("ANOVA test:", sep="\n")
anovaTest <- aov(values ~ class, data = rates_t)
print(TukeyHSD(x = anovaTest, 'class', conf.level = 0.95)$class)
cat("\n")
cat(paste0("Saving output:\n\t- ", outputPDF), sep="\n")

pdf(file=outputPDF, width=8, height=6)
                print(rates_p)
        invisible(dev.off())
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


