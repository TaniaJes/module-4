#!/bin/bash 
#SBATCH --job-name=Kallisto # Job name 
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=tania.gonzalez-robles@nyulangone.org # Where to send mail 
#SBATCH --ntasks=6 # Run on a single CPU 
#SBATCH --mem=40gb # Job memory request 
#SBATCH --time=10:00:00 # Time limit hrs:min:sec 
#SBATCH --output=/gpfs/home/tgr248/jobreports/slurm_%A_%a.log # Standard output and error log 
#SBATCH -p cpu_short # Specifies location to submit job
#SBATCH --array=609-612,615,616


# LOAD MODULES
module load subread/1.6.3     # contains featureCounts command
module load kallisto/0.44.0


# ENVIRONMENT SETUP
mkdir -p /gpfs/scratch/tgr248/bioinformatics # -p flag will create a folder if necessary
cd /gpfs/scratch/tgr248/bioinformatics
mkdir -p ./module-3 ./module-3/data ./module-4 ./module-4/data


DATA_ID="7049" # beginning of SRR ID, without the array section

DIRECTORY=/gpfs/scratch/tgr248/bioinformatics/module-3/data
OUTPUT_DIR=/gpfs/scratch/tgr248/bioinformatics/module-4/data
ANNOTATION_FILE=/gpfs/data/courses/bminga3004/Practicum5/genes.gtf
TRANSCRIPT_REF=/gpfs/data/courses/bmscga2604/gencode.v34.pc_transcripts.fa # assuming human transcript
SPECIE=HomoSapiens


# JOB DESCRIPTION
hostname
date
pwd
echo "SRR$DATA_ID$SLURM_ARRAY_TASK_ID" # to track which is current file 
mkdir -p $OUTPUT_DIR/featureCounts_output/SRR${DATA_ID}${SLURM_ARRAY_TASK_ID}

# GENE COUNTS
featureCounts -s2 -p -T6 -a $ANNOTATION_FILE -o $OUTPUT_DIR/featureCounts_output/SRR${DATA_ID}${SLURM_ARRAY_TASK_ID}/featureCounts.output.txt $DIRECTORY/samtools/SRR${DATA_ID}${SLURM_ARRAY_TASK_ID}_trimmed_bt2_mapped.sorted.bam

#-s2        Perform strand-specific read counting. Usually needed when treating mRNA U degradase.
#           	Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
# -p        If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads.
# -T        Number of the threads. 1 by default.
# -a        Name of an annotation file. GTF/GFF
# -o        Name of the output file including read counts.
# <input>    A list of SAM or BAM format files. Sorted. 



# TRANSCRIPT COUNTS FROM PSEUDOALIGNMENT (KALLISTO)

### Indexing (only need to do this once)
kallisto index -i $OUTPUT_DIR/${SPECIE}.transcript.index $TRANSCRIPT_REF

### transcrip quantification
kallisto quant -i $OUTPUT_DIR/${SPECIE}.transcript.index -o $OUTPUT_DIR/SRR${DATA_ID}${SLURM_ARRAY_TASK_ID} -b 100 --bias --fr-stranded $DIRECTORY/SRR${DATA_ID}${SLURM_ARRAY_TASK_ID}_1.fastq.gz $DIRECTORY/SRR${DATA_ID}${SLURM_ARRAY_TASK_ID}_2.fastq.gz
# Note: According to class slides, I should feed untrimmed fastq files to pseudoaliners like Kallisto.
# default paired-end mode

# -b (100) 	bootstrap
# --bias	performs sequence-based bias correction
# --ft-stranded	strandness specification (options: --ft-stranded, --tf-stranded, --ft-unstranded)


### END OF CODE ###
