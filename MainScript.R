if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("janitor")
install.packages("pheatmap")
install.packages("ggrepel")

BiocManager::install("GEOquery")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")

library(GEOquery)
library(conflicted)
library(dplyr)
library(tibble)
library(janitor)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(WGCNA)
library(biomaRt)

#gse <- getGEO("GSE64810")
#diagnosis <- gse[[1]]$`diagnosis:ch1`
#disease_label <- ifelse(diagnosis == "Neurologically normal", FALSE, TRUE)

expressions <- read.table("GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt.gz", header=TRUE, sep = "\t", comment.char = "!")
counts <- read.table("GSE64810_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")

padj_threshold <- 0.01
log2FoldChange_threshold = 1

expressions$level <- "Not DE"
expressions$level[expressions$log2FoldChange >= log2FoldChange_threshold & expressions$padj < padj_threshold] <- "Over"
expressions$level[expressions$log2FoldChange <= -log2FoldChange_threshold & expressions$padj < padj_threshold] <- "Under"

overexpressed <- expressions[expressions$level == "Over",]
underexpressed <- expressions[expressions$level == "Under",]

# order sets of over- and underexpressed genes in a descending fashion
overexpressed <- overexpressed[order(-overexpressed$log2FoldChange),]
underexpressed <- underexpressed[order(underexpressed$log2FoldChange),]

# pick the topmost over- and underexpressed genes to depict them in the volcano plot following
selection <- 10
top_over <- head(overexpressed, selection)
top_under <- head(underexpressed, selection)
top_degs <- rbind(top_over, top_under)

colours <- c("Over" = "#de28c3", "Not DE" = "#b7de28", "Under" = "#369bba")
sizes <- c("Over" = 1, "Not DE" = 0.7, "Under" = 1)
alphas <- c("Over" = 0.5, "Not DE" = 0.1, "Under" = 0.5)
ggplot(data=DataFrame(expressions), aes(x=log2FoldChange, y=-log10(padj), fill=level, size=level, alpha=level, label=symbol)) + 
  geom_point(shape=21,colour="black") + 
  scale_fill_manual(values=colours) + 
  scale_size_manual(values=sizes) + 
  scale_alpha_manual(values=alphas) +
  geom_label_repel(data = top_degs, aes(label = symbol), force = 2, nudge_y = 1, size = 2)

mart <- useMart("ENSEMBL_MART_ENSEMBL");
mart <- useDataset('hsapiens_gene_ensembl', mart)
hgnc <- getBM(mart = mart, values = expressions[expressions$padj <= 0.05,]$symbol, attributes=c("entrezgene_id"), filter="external_gene_name")

c <- dplyr::filter(counts, GeneID %in% hgnc$entrezgene_id)
temp <- c[,-1];
rownames(temp) <- c$GeneID;
c <- as.data.frame(temp);

powers = c(1:20)
gsg <- goodSamplesGenes(c, verbose=3)
if (!gsg$allOK) {
  counts <- counts[gsg$goodGenes, gsg$goodSamples]
}
sft = pickSoftThreshold(t(c), powerVector=powers, verbose=5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type="n", xlab="Soft threshold (power)", ylab="Scale Free Topology model fit")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.8, col="red")

softPower = 11
conflicts_prefer(WGCNA::cor)
adjacency = adjacency(t(c), power = softPower)

# Convert adjacency to Topological Overlap Matrix (TOM) and calculate the dissimilarity
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# Perform hierarchical clustering based on TOM dissimilarity
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based Dissimilarity", labels=FALSE, hang=0.04)

# Dynamically cut tree to identify modules
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods)

# Convert numeric module labels to colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Calculate module eigengenes
MEList = moduleEigengenes(t(c), colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate module membership (kME) for each gene
moduleMembership = cor(t(c), MEs, use = "pairwise.complete.obs")
moduleColors = dynamicColors

# Focus on a specific module (e.g., "blue")
module = "blue"
modGenes = (moduleColors == module)

# Get gene names and their connectivity in the module
hubGenes = names(which.max(moduleMembership[modGenes, "MEblue"]))
cat("Hub gene for module", module, ":", hubGenes, "\n")
