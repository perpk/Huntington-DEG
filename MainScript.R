if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("janitor")

BiocManager::install("GEOquery")
BiocManager::install("DESeq2")

library(GEOquery)
library(conflicted)
library(dplyr)
library(tibble)
library(janitor)
library(DESeq2)
library(ggplot2)

gse <- getGEO("GSE64810")
diagnosis <- gse[[1]]$`diagnosis:ch1`
disease_label <- ifelse(diagnosis == "Neurologically normal", FALSE, TRUE)

expressions <- read.table("GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz", header=TRUE, sep = "\t", comment.char = "!")
# expressions <- data.frame(expressions[,-1], row.names=expressions[,1])
# expressions_t <- setNames(as.data.frame(t(expressions[-1])), expressions[[1]])
# expressions_t <- expressions_t %>% mutate(status=disease_label)

# expressions$diffxpressd <- expressions$padj < 0.01
padj_threshold <- 0.01
log2FoldChange_threshold = 1

expressions$level <- "NA"
expressions$level[expressions$log2FoldChange >= log2FoldChange_threshold & expressions$padj < padj_threshold] <- "OVER"
expressions$level[expressions$log2FoldChange <= -log2FoldChange_threshold & expressions$padj < padj_threshold] <- "UNDER"

colours <- c("OVER" = "#de28c3", "NA" = "#b7de28", "UNDER" = "#369bba")
sizes <- c("OVER" = 1, "NA" = 0.7, "UNDER" = 1)
alphas <- c("OVER" = 0.5, "NA" = 0.1, "UNDER" = 0.5)
ggplot(data=DataFrame(expressions), aes(x=log2FoldChange, y=-log10(padj), fill=level, size=level, alpha=level)) + geom_point(shape=21,colour="black") + scale_fill_manual(values=colours) + scale_size_manual(values=sizes) + scale_alpha_manual(values=alphas)

