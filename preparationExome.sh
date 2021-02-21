#!/bin/bash
# Exome analysis for detecting somatic variants
# script 1: preparationExome.sh   
#--------------------------------------------------------------------------
#This script does the following:
#installs the necessary tools and their dependencies
#creates working direcrories
#downloads and decompresses the input fastq files
#downloads the reference genome and performs its indexing with BWA
#downloads and decompresses the GTF annotation;

#When this script is run, please, proceed with script 2: pipelineExome.sh
#--------------------------------------------------------------------------

###For commenting out pieces of the code
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

###Introducing variables for paths so we can change them later and reuse the script
path_to_project="/home/ubuntu/data/mydatalocal/TP2_exome/"    
path_to_raw_fastq="$path_to_project""raw_fastq/"
path_to_index="$path_to_project""bwa_index/"
path_to_gtf="$path_to_project""gtf/"
genome="chr16.fa"

###Installing the tools
echo -e "Installing the necessary tools\n"
conda activate
sudo apt install fastqc   #for quality control
conda install trimmomatic  #for trimming the reads
sudo apt install bwa  #for genome indexing and read mapping with BWA
sudo apt install samtools #for SRA processing
sudo apt install varscan #for variants calling
sudo apt install default-jre #dependency package pointing to the Java runtime
sudo apt install bedtools  #for annotation part 
echo -e "\nAll installations completed!\n"

###Creating directories
echo "Creating project directories"
mkdir "$path_to_project"
mkdir "$path_to_raw_fastq"
mkdir "$path_to_index"
mkdir "$path_to_gtf"
cd "$path_to_project"
echo -e "Project directories created\n"

###Download and decompress fastq files
echo -e "Downloading raw fastq files\n"
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt 
tar -zxvf patient7.tar.gz --directory "$path_to_raw_fastq" #2>/dev/null 
echo -e "\nFiles downloaded to " "$path_to_raw_fastq"patient7.exome/

path_to_raw_fastq="$path_to_raw_fastq""patient7.exome/"

###Downloading and indexing the genome 
echo -e "\nDownloading reference genome"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz -q
gunzip chr16.fa.gz
mv ./chr16.fa "$path_to_index"
echo "Genome downloaded to " "$path_to_index"
echo -e "\nCreating BWA index\n"
bwa index -a bwtsw "$path_to_index"chr16.fa
echo -e "\nIndex created!\n"

###Downloading GTF annotation for annotating VCF
echo -e "\nDownloading GTF annotation"
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz -q
gunzip gencode.v24lift37.basic.annotation.gtf.gz
mv ./gencode.v24lift37.basic.annotation.gtf "$path_to_gtf"
echo -e "\nGTF downloaded to " "$path_to_gtf"

echo "Everything is prepared. You are ready to run script2: pipelineExome.sh"
