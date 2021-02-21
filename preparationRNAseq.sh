#!/bin/bash
# RNA-seq analysis 
# script 1: preparationRNAseq.sh   
#--------------------------------------------------------------------------
#This script does the following:
#installs the necessary tools and their dependencies
#creates working direcrories
#downloads and decompresses the input fastq files
#downloads and decompresses the reference genome and the GTF annotation
#creates STAR genome index

#When this script is run, please, proceed with script 2: pipelineRNAseq.sh
#--------------------------------------------------------------------------


###For commenting out pieces of the code
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

###Introducing variables for paths so we can change them later and reuse the script
path_to_project="/home/ubuntu/data/mydatalocal/TP1_RNAseq/"
path_to_raw_fastq="$path_to_project""raw_fastq/"
path_to_ref_genome="$path_to_project""ref_genome/"
path_to_index="$path_to_project""STARindex/"
annotation="annotation.gtf"
genome="chr18.fa"

###Installing the tools
echo "Installing the necessary tools"
conda activate
sudo apt install fastqc #for quality control
conda install trimmomatic #for trimming the reads
sudo apt install rna-star #for genome indexing and read mapping with STAR
sudo apt install subread #for counting reads via featureCounts
sudo apt install samtools #for indexing bam
echo "All installations completed!"

###Creating directories
echo "Creating project directories"
mkdir "$path_to_project"
mkdir "$path_to_ref_genome"
mkdir "$path_to_raw_fastq"
mkdir "$path_to_index"
cd "$path_to_project"
echo "Project directories cerated"

###Download and decompress fastq files
echo "Downloading raw fastq files"
wget  http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz -q
tar -zxvf TPrnaseq.tar.gz --directory "$path_to_raw_fastq"
echo "Files downloaded to " "$path_to_raw_fastq"

###Download and decompress the genome and annotation
echo "Downloading reference genome"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz -q
gunzip chr18.fa.gz
mv ./chr18.fa "$path_to_ref_genome$genome"
echo "Genome downloaded to " "$path_to_ref_genome"

echo "Downloading annotation"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz -q
gunzip gencode.v24lift37.basic.annotation.gtf.gz
mv ./gencode.v24lift37.basic.annotation.gtf "$path_to_ref_genome$annotation"
echo "Annotation downloaded to " "$path_to_ref_genome"

###Creating STAR genome index
echo "Starting to build STAR genome index"

STAR --runMode genomeGenerate --runThreadN 10 \
  --genomeDir "$path_to_index" \
  --genomeFastaFiles "$path_to_ref_genome$genome" \
  --sjdbGTFfile "$path_to_ref_genome$annotation"

echo -e "\nSTAR index saved to " "$path_to_index"

echo "Everything is prepared. You are ready to run script2: pipelineRNAseq.sh"
