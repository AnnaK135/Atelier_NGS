#!/bin/bash
#!/bin/bash
# RNA-seq analysis 
# script 3: quality_controlRNAseq.sh
#--------------------------------------------------------------------------
# This script performs an optional quality control step for fastq files.
# Should be run no earlier than script1: preparationRNAseq.sh
#--------------------------------------------------------------------------


###for commenting out pieces of the code
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

path_to_project="/home/ubuntu/data/mydatalocal/TP1_RNAseq/"
path_to_raw_fastq="$path_to_project""raw_fastq/"
path_to_trimmed_fastq="$path_to_project""trimmed_fastq"

### Assessing raw reads
echo "Quality assessment for the raw reads:"
echo "Checking for the necessary folders"
if [ ! -d "$path_to_raw_fastq" ];then
  echo "Necessary folders do not exist! Run the script preparationTP.sh first"
  exit 1
fi
echo "Folders are OK"

cd "$path_to_project"
mkdir fastqc
fastqc "$path_to_raw_fastq"*.fastq --outdir "$path_to_project"fastqc
echo "Quality check done! Report files are saved to" "$path_to_project""fastqc"

### Assessing trimmed reads, if trimming has already been done
if [ ! -d "$path_to_trimmed_fastq" ];then
  echo "Trimming has not been done. Nothing more to assess. Exiting"
  exit 0
fi
echo "Quality assessment for the trimmed reads:"
fastqc "$path_to_trimmed_fastq"*.fastq --outdir "$path_to_project"fastqc
echo "Quality check for the trimmed reads done! Report files are saved to" "$path_to_project""fastqc"
