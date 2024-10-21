BiocManager::install("Biobase")
BiocManager::install("GO.db")
BiocManager::install("impute")
install.packages("WGCNA")

library(WGCNA)

counts <- read.table("GSE64810_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")
# gse <- getGEO("GSE64810")
# diagnosis <- gse[[1]]$`diagnosis:ch1`
# disease_label <- ifelse(diagnosis == "Neurologically normal", FALSE, TRUE)

temp <- counts[,-1];
rownames(temp) <- counts$GeneID;
counts <- as.data.frame(temp);

powers = c(1:20)
disableWGCNAThreads()
gsg <- goodSamplesGenes(counts, verbose=3)
if (!gsg$allOK) {
  counts <- counts[gsg$goodGenes, gsg$goodSamples]
}
gc()
install.packages("bigmemory")
install.packages("biganalytics")
library(bigmemory)
library(biganalytics)

# cpu/time costly operation.- runs over >30k genes from the dataset.
counts_bm <- as.big.matrix(t(counts), backingfile="big_counts.bin", descriptorfile="big_counts.desc")

sft = pickSoftThreshold(t(counts), powerVector=powers, verbose=5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type="n", xlab="Soft threshold (power)", ylab="Scale Free Topology model fit")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")

softPower = 6
conflicts_prefer(WGCNA::cor)
enableWGCNAThreads(nThreads = 4)
adjacency = adjacency(t(counts), power = softPower)

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
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate module membership (kME) for each gene
moduleMembership = cor(datExpr, MEs, use = "pairwise.complete.obs")
moduleColors = dynamicColors

# Focus on a specific module (e.g., "blue")
module = "blue"
modGenes = (moduleColors == module)

# Get gene names and their connectivity in the module
hubGenes = names(which.max(moduleMembership[modGenes, module]))
cat("Hub gene for module", module, ":", hubGenes, "\n")
