# Install WGCNA if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("WGCNA")

# Load the WGCNA package
library(WGCNA)

# Simulate random gene expression data for demonstration
set.seed(123)
datExpr <- matrix(rnorm(1000*20), nrow=1000, ncol=20) # 1000 genes, 20 samples
rownames(datExpr) <- paste("Gene", 1:1000, sep="")
colnames(datExpr) <- paste("Sample", 1:20, sep="")

# Check for missing data
gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodGenes, gsg$goodSamples]
}

# Choose the soft-thresholding power
powers = c(1:20)
sft = pickSoftThreshold(t(datExpr), powerVector = powers, verbose = 5)

# Plot the scale-free topology fit index as a function of the power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type="n", 
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit", 
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")

# Set soft-thresholding power based on previous step (e.g., power = 6)
softPower = 6
conflicts_prefer(WGCNA::cor)
adjacency = adjacency(t(datExpr), power = softPower)

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
