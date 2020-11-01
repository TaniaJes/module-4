### LOAD REQUIRED LIBRARIES
library(biomaRt)
library(tximport)
library(rhdf5)

### SET WORKING DIRECTORY ### You will need to edit this and direct it your downloaded kallisto folder
setwd("/Users/taniagonzalezrobles/NYU-classes/bioinformatics/module-4/data/kallisto_output")

### IMPORT ENSEMBl ANNOTATIONS FOR HUMAN GENOME & GENERATE TWO COLUMN FILE LINKING TRANSCRIPT AND GENE IDS
mart <- biomaRt::useMart(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype", "refseq_mrna", "refseq_ncrna"), mart = mart)
t2g$target_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep=".") # append version number to the transcript ID
t2g[,c("ensembl_transcript_id","transcript_version")] <- list(NULL) # delete the ensembl transcript ID and transcript version columns
t2g <- dplyr::rename( t2g, gene_symbol = external_gene_name, full_name = description, biotype = transcript_biotype ) # to add gene symbols
t2g<-t2g[,c(ncol(t2g),1:(ncol(t2g)-1))]

### GENERATE SECOND OBJECT CONTAINING ONLY PROTEIN CODING GENES
gb <- getBM(attributes=c("ensembl_gene_id","gene_biotype"), mart=mart)
gb_coding<-subset(gb, gb$gene_biotype=="protein_coding")
genes<-gb_coding$ensembl_gene_id

### USE TXIMPORT TO SUMMARIZE TRANSCRIPT COUNTS INTO GENE COUNTS

## For multiple samples, each named as a folder in the kallisto directory (can be abundance.h5 or abundance.tsv file)
accessions <- list.dirs(full.names=FALSE)[-1]
kallisto.dir<-paste0(accessions)
kallisto.files<-file.path(".",kallisto.dir,"abundance.h5") #can also be abundance.tsv
names(kallisto.files)<- accessions
tx.kallisto <- tximport(kallisto.files, type = "kallisto", tx2gene = t2g, countsFromAbundance ="no", ignoreAfterBar=T) # ignoreAfterBar=T will omit the | separator in file ID.

### GENERATE TWO COLUMN OUTPUT FORMAT
counts<-as.data.frame(tx.kallisto$counts[row.names(tx.kallisto$counts) %in% genes, ])
ids<-rownames(counts)

### ROUND VALUES (DESEQ2 DOES NOT LIKE FRACTIONS), AND WRITE TO OUTPUT FILE
write.table(round(counts),paste("kallisto-gene-output",".txt",sep=""), row.names=ids, quote=F, col.names=T, sep="\t")


### END OF CODE ###

