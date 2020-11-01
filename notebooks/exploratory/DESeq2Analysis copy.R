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
library("EnhancedVolcano")
library("RColorBrewer")
library("biomaRt")

### SET WORKING DIRECTORY
# note: this directory should be populated with the raw counts file
setwd("~/NYU-classes/bioinformatics/module-4/data")


### Import count table and details on experimental design
# NB: Make sure column names in the sample(table) file and counts file are exactly the same and in the same order
samples <- read.table("~/NYU-classes/bioinformatics/module-4/data/sample_info.txt", header = TRUE)
featCounts <- read.table("~/NYU-classes/bioinformatics/module-4/data/featureCounts_output/geneCounts-output.txt", header = T, row.names = 1)
featCounts <- featCounts[, rownames(samples)] # column reordering to match samples order
Dataset <- DESeqDataSetFromMatrix(countData = featCounts, colData = samples, design = ~batch + condition)
# Note: Always end with conditions for 'design' variable


### PRELIMINARY ANALYSES ###
# The first steps in your analysis should focus on better understanding the relationship of the datasets being studied. 
# This can be simply achieved by generating a PCA plot showing the relationship of your samples.
# First we transform our raw count data using a variance stabilizing transformation (VST) that roughly mirrors how DeSeq2 models the data.
vsd1 <- varianceStabilizingTransformation(Dataset, blind=FALSE)

# Then we plot a PCA, grouping and coloring our datasets according to batch
plotPCA(vsd1, "condition")
### note that you can attach additional information based on the column headers in your sample table
plotPCA(vsd1, c("condition","batch"))

# we can also attempt to replicate the batch effect correction performed by DeSeq2 using the limma::removeBatchEffect function
vsd2 <- varianceStabilizingTransformation(Dataset, blind=FALSE)
assay(vsd2) <- limma::removeBatchEffect(assay(vsd2), vsd2$batch)
plotPCA(vsd2, "condition")