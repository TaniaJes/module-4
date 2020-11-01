
### LOAD REQUIRED LIBRARIES
library(biomaRt)
library(tximport)
library(rhdf5)

### SET WORKING DIRECTORY ### 
setwd("/Users/taniagonzalezrobles/NYU-classes/bioinformatics/module-4/data/featureCounts_output")

### IMPORT ENSEMBL ANNOTATIONS FOR HUMAN GENOME & GENERATE TWO COLUMN FILE LINKING TRANSCRIPT AND GENE IDS
mart <- biomaRt::useMart(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype", "refseq_mrna", "refseq_ncrna"), mart = mart)
t2g$target_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep=".") # append version number to the transcript ID
t2g[,c("ensembl_transcript_id","transcript_version")] <- list(NULL) # list(NULL) deletes the ensembl transcript ID and transcript version columns
t2g <- dplyr::rename( t2g, gene_symbol = external_gene_name, full_name = description, biotype = transcript_biotype )
t2g<-t2g[,c(ncol(t2g),1:(ncol(t2g)-1))] #to move columns order

### GENERATE ADDITIONAL OBJECT CONTAINING ONLY PROTEIN CODING GENES
gb <- getBM(attributes=c("ensembl_gene_id","gene_biotype"), mart=mart)
gb_coding <- subset(gb, gb$gene_biotype=="protein_coding")
genes <- gb_coding$ensembl_gene_id # create new elements with the list of genes_id for protein coding genes only

### featureCounts files import (.txt file format)

# Import of multiple datasets after counting with featureCount through BigPurple
bt2.files <- list.files(full.names=F)
bt2.path <- file.path(".", bt2.files, "featureCounts.output.txt")

### GENERATE TWO COLUMN OUTPUT FORMAT
# set counts matrix
counts <- c()
for( i in seq_along(bt2.path) ){
  x <- read.table(file=bt2.path[i], sep="\t", header=T)
  counts <- cbind(counts, x[,7]) # V7 is the counts column
}

ids <- gsub("\\..*","", x[,1]) # V1 is the Geneid, gsub: to remove extra numbers after the "." in ENSBL IDs
rownames(counts) <- ids # setting row names
colnames(counts) <- bt2.files # setting column names

# cut -f1,7 ./SRR.../featureCounts.output.txt > geneCounts.txt # can use something like this to extract the gene counts from terminal

counts <- as.data.frame(counts[row.names(counts) %in% genes, ]) # filtering only protein coding genes
ids <- rownames(counts)

### ROUND VALUES (DESEQ2 DOES NOT LIKE FRACTIONS), AND WRITE TO OUTPUT FILE
write.table(round(counts),paste("geneCounts-output",".txt",sep=""), row.names=ids, quote=F, col.names=T, sep="\t")

### END OF CODE ###

