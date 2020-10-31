##########################################################################################
# This script performs differential gene expression analysis with DESeq2 #
##########################################################################################

### LOAD REQUIRED LIBRARIES
library("DESeq2")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("ggplot2")
library("ashr") # library for type = "ashr" object used to generate DEseq results object
library("gridExtra")
library("plotly")
library("pheatmap")
library("EnhancedVolcano") # pretty volcano plots
library("RColorBrewer")
library("biomaRt")

### SET WORKING DIRECTORY
# note: this directory should be populated with the raw counts file
setwd("~/NYU-classes/bioinformatics/module-4/data")


### Import count table and details on experimental design
# NB: Make sure column names in the sample(table) file and counts file are exactly the same and in the same order
samples <- read.table("./sample_info.txt", header = TRUE)
kallistoCounts <- read.table("./kallisto_output/kallisto-gene-output.txt", header = T, row.names = 1)
kallistoCounts <- kallistoCounts[, rownames(samples)] # column reordering to match samples order
Dataset <- DESeqDataSetFromMatrix(countData = kallistoCounts, colData = samples, design = ~batch + condition)
# Note: Always end with conditions for 'design' variable


### PRELIMINARY ANALYSES ###
# First we transform our raw count data using a variance stabilizing transformation (VST) that roughly mirrors how DeSeq2 models the data.
vst1 <- varianceStabilizingTransformation(Dataset, blind = FALSE)
plotPCA(vst1, "condition") # Then we plot a PCA, grouping and coloring our datasets according to batch
plotPCA(vst1, c("condition","batch")) 
# note: you can attach additional information based on the column headers in your sample table

# We can also attempt to replicate the batch effect correction performed by DeSeq2 using the `limma::removeBatchEffect` function.
vst2 <- varianceStabilizingTransformation(Dataset, blind = FALSE)
assay(vst2) <- limma::removeBatchEffect(assay(vst2), vst2$batch)
plotPCA(vst2, "condition")
plotPCA(vst2, c("condition","batch"))

# uncorrected
sampleDists <- dist( t( assay(vst1) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main = "Samples without batch correlation")

# corrected
sampleDistsCorr <- dist( t( assay(vst2) ) )
sampleDistsCorr
sampleDistCorrMatrix <- as.matrix( sampleDistsCorr )
colnames(sampleDistCorrMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
pheatmap(sampleDistCorrMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, main = "Samples with batch correlation")

# At this stage, you should have a good sense of how your samples cluster and the effect of batch correction (if used)
# In simple RNA-Seq situations (control vs treatment, 3-5 bioreps each), this is all that should be required. 
# For more complex situations, you will need to dive deep into working of DeSeq2
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


### BASIC DGE ANALYSIS USING DESEQ2 ###

# Run DESEQ and generate a simple plot showing the distribution of regulated and unregulated genes
DatasetProcessed <- DESeq(Dataset) # runs DESEQ
plotMA(DatasetProcessed, main="Heteroskedasticity assessment", ylim=c(-8,15)) # red dots have p val < 0.05

# Next we perform a contrast analysis to produce a list of differentially regulated genes between our two conditions
# First we set CTRL dataset as baseline
Dataset$condition <- relevel(Dataset$condition, "siCtrlBH5")

# Next we create our results object while performing shrinkage of effect size 
# (this reduces the impact of apparent gross changes in low expressed genes)
res1 <- lfcShrink(DatasetProcessed, contrast=c("condition","siCtrlBH5_dsDNA","siCtrlBH5"), type = "ashr")
plotMA(res1, main="Restored heteroskedasticity", ylim=c(-8,15))

# Data summary
summary(results(DatasetProcessed, alpha = 0.01, lfcThreshold = 1))
# out of 15324 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 1.00 (up)    : 115, 0.75%
# LFC < -1.00 (down) : 16, 0.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 10237, 67%
# (mean count < 7)

# Here we modify our output data to include two additional columns that contain the baseMeans (a proxy for counts)
# This is useful for downstream filtering of lowly expressed genes
baseMean.siCtrlBH5 <-  rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "siCtrlBH5"])
baseMean.siCtrlBH5_dsDNA <-  rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "siCtrlBH5_dsDNA"])

res1 <-  cbind(as.data.frame(res1), baseMean.siCtrlBH5, baseMean.siCtrlBH5_dsDNA)

# Here we add two further columns, the gene symbol (common name) and entrez ID - both of which may be useful downstream
res1$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="SYMBOL", keytype="ENSEMBL", multiVals="first") # MAPS GENE IDs
res1$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Finally we write the complete results object to an outfile
write.csv(res1, "kallisto_DESeq2.csv", row.names=TRUE)


### PLOTTING AND DATA VISUALIZATION ###

## Volcano plot
EnhancedVolcano(toptable = res1, x = "log2FoldChange", y = "pvalue",
                lab = res1$symbol, pCutoff = 10e-6, FCcutoff = 1, xlim = c(-5, 5), ylim = c(-20, 100),
                title = "Kallisto DESeq2 results", subtitle = "Differential Expression Analysis", caption = NULL,
                col = c("grey30", "forestgreen", "purple", "dark orange"))


## Plot heatmap - top20 variable genes
rld <- rlog(DatasetProcessed) # apply "regularized log" transformation to DESeqDataSet
top20Genes <- order(rowMeans(assay(rld)), decreasing = T)[1:20] # Extract (assay) matrix for top 20 variable genes
# bottom20Genes <- 
df <- as.data.frame(colData(DatasetProcessed)) # create annotation matrix: c("batch", "condition", "sizeFactor")
mat <- assay(rld)[top20Genes, ] # generate new matrix w/top20 variable genes
mat<-mat-rowMeans(mat) # plot distances from the mean, makes heatmap clearer (why is this important?)
mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast") 
# to circunvent server overloading error 
# assuming human, to change geneID to symbols. 
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1] # to change geneIDs to symbol in "mat" matrix.
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, main = "Kallisto - Top 20 Variable Genes")

# Plot heatmap - sorting by padj and log2fc
# PICK ALL GENES WITH pADJ < 0.01 AND THEN SUBSET FOR THOSE WITH Log2FC > 1 THEN PICK TOP 25 HITS
subset.top <- head((subset(res1, res1$padj < 0.01)), n=25)  
# subset.bottom <- tail((subset(res1, res1$padj < 0.01)), n=25)
sigGenes <- rownames(subset.top)
rows <- match(sigGenes, row.names(rld))
mat <- assay(rld)[rows,] # to extract matrix 
mat <- mat-rowMeans(mat) ## plot distances from the mean, makes heatmap clearer
mart <- useMart("ensembl","hsapiens_gene_ensembl") ## assuming human
gns <- getBM(c("hgnc_symbol","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1]
df <- as.data.frame(colData(DatasetProcessed))
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, main = "Kallisto - Top 25 most significant genes") 


### END OF CODE ###

