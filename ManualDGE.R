BiocManager::install("biomaRt")

require(biomaRt)

counts <- read.table("GSE64810_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")
gse <- getGEO("GSE6x4810")
diagnosis <- gse[[1]]$`diagnosis:ch1`
disease_label <- ifelse(diagnosis == "Neurologically normal", FALSE, TRUE)

mart <- useMart("ENSEMBL_MART_ENSEMBL");
mart <- useDataset('hsapiens_gene_ensembl', mart)
hgnc <- getBM(mart = mart, values = counts$GeneID, attributes = c('hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id'), filter="entrezgene_id")

temp <- counts[,-1];
rownames(temp) <- counts$GeneID;
counts <- as.data.frame(temp);


