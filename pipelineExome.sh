#!/bin/bash
# Exome analysis for detecting somatic variants
# script 2: pipelineExome.sh   
#--------------------------------------------------------------------------
# This script does the following:
# Trims fastq files to eliminate low quality bases at read extremities
# Maps the trimmed reads to the indexed reference genome -> .sam
# Collects mapping statistics (.stats)
# Converts SAM files to sorted BAM files an then to PileUp files
# Performs BAM indexing (for visualisation in IGV)
# Calls somatic variants from Normal and Tumor PileUp files -> .vcf
# Extracts somatic mutations from a VCF file and performs basic annotation 

#Before running this script, please, run the script 1: preparationExome.sh
#--------------------------------------------------------------------------


###For commenting out pieces of the code
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

###Introducing variables for paths so we can change them later and reuse the script
path_to_project="/home/ubuntu/data/mydatalocal/TP2_exome/"  
path_to_raw_fastq="$path_to_project""raw_fastq/"
path_to_raw_fastq="$path_to_raw_fastq""patient7.exome/"
path_to_index="$path_to_project""bwa_index/"
path_to_gtf="$path_to_project""gtf/"
genome="chr16.fa"
gtf="gencode.v24lift37.basic.annotation.gtf"
normal="TCRBOA7-N-WEX-chr16"
tumor="TCRBOA7-T-WEX-chr16"

###Sanity checks
echo "Checking for the necessary folders"
if [ ! -d "$path_to_project" ];then
  echo "Necessary folders do not exist! Please, run the script preparationExome.sh first"
  exit 1
fi

if [ ! -d "$path_to_raw_fastq" ];then
  echo "Necessary folders do not exist! Please, run the script preparationExome.sh first"
  exit 1
fi

if [ ! -d "$path_to_gtf" ];then
  echo "GTF annotation file does not exist! Please, run the script preparationExome.sh first"
  exit 1
fi

if [ ! -d "$path_to_index" ];then
  echo "BWA index does not exist. Please, run preparationExome.sh"
  exit 1
fi
echo -e "Folders are OK\n"

echo "Checking for the index files"
if [ ! "$path_to_index""$genome".bwt ];then
  echo "BWA index does not exist or is incomplete. Please, run preparationExome.sh"
  exit 1
fi

if [ ! "$path_to_index""$genome".amb ];then
  echo "BWA index does not exist or is incomplete. Please, run preparationExome.sh"
  exit 1
fi

if [ ! "$path_to_index""$genome".ann ];then
  echo "BWA index does not exist or is incomplete. Please, run preparationExome.sh"
  exit 1
fi

if [ ! "$path_to_index""$genome".pac ];then
  echo "BWA index does not exist or is incomplete. Please, run preparationExome.sh"
  exit 1
fi

if [ ! "$path_to_index""$genome".sa ];then
  echo "BWA index does not exist or is incomplete. Please, run preparationExome.sh"
  exit 1
fi
echo -e "The index exists. Proceeding...\n"


###Trimming low quality bases at read extermities
echo "Trimming fastq files"
ls "$path_to_raw_fastq" | sed -e 's/_r[1,2]F.fastq.gz$//' | uniq > "$path_to_project"fastq_prefixes
mkdir "$path_to_project"trimmed_fastq

for i in $(cat "$path_to_project"fastq_prefixes); do
  trimmomatic PE "$path_to_raw_fastq"${i}_r1F.fastq.gz "$path_to_raw_fastq"${i}_r2F.fastq.gz -baseout "$path_to_project"trimmed_fastq/${i}.trimmed.fastq LEADING:20 TRAILING:20 MINLEN:50
done

rm "$path_to_project"trimmed_fastq/*U.fastq #removing the files with trimmed low quality reads
path_to_trimmed_fastq="$path_to_project""trimmed_fastq/"

if [ ! -d "$path_to_trimmed_fastq" ];then
	echo "Trimming failed. Exiting"
	exit 1
fi

###Mapping the reads with BWA
echo -e "Started BWA mapping...\n"
mkdir "$path_to_project"sam
path_to_sam="$path_to_project""sam/"
for i in $(cat "$path_to_project"fastq_prefixes); do
    bwa mem -M -t 2 -A 2 -E 1 "$path_to_index""$genome" "$path_to_trimmed_fastq"${i}.trimmed_1P.fastq "$path_to_trimmed_fastq"${i}.trimmed_2P.fastq > \
    "$path_to_sam"${i}.sam
done
echo -e "Mapping done! SAM files are saved to" "$path_to_sam"


###Processing SAM files
echo -e "\nStarted processing SAM files"
path_to_bam="$path_to_project""bam/"
path_to_mpileup="$path_to_project""mpileup/"
mkdir "$path_to_project"bam
mkdir "$path_to_project"mpileup

for i in $(cat "$path_to_project"fastq_prefixes); do
    echo -e "\nStarted processing sample" "${i}"
    echo -e "Converting SAM to BAM\n"
    samtools view -b "$path_to_sam"${i}.sam -o "$path_to_bam"${i}.bam
    echo -e "Sorting BAM\n"
    samtools sort "$path_to_bam"${i}.bam -o "$path_to_bam"${i}.bam
    echo -e "Indexing BAM\n"
    samtools index "$path_to_bam"${i}.bam
    echo -e "Collecting mapping statistics"
    samtools flagstat "$path_to_bam"${i}.bam > "$path_to_bam"${i}.stats
    echo -e "Converting BAM to Mpileup\n"
    samtools mpileup -B -A -f "$path_to_index""$genome" "$path_to_bam"${i}.bam > "$path_to_mpileup"${i}.pileup
done
echo -e "Processing done! Stats, BAM and BAI files are saved to" "$path_to_bam"
echo -e "PileUp files saved to" "$path_to_mpileup"


###Calling somatic variants with Varscan
echo -e "\nCalling somatic variants"
path_to_vcf="$path_to_project""vcf/"
mkdir "$path_to_project"vcf

varscan somatic "$path_to_mpileup""$normal".pileup "$path_to_mpileup""$tumor".pileup "$path_to_vcf"variants --variants --p-value 0.001 --min-avg-qual 15 --output-vcf 1
echo -e "\nDone! VCF is saved to" "$path_to_vcf"
  
###VCF basic annotation
echo -e "\nExtracting and annotating somatic mutations"
path_to_bed="$path_to_project""bed/"
path_to_annotated_somatic_mutations="$path_to_project""annotated_somatic_mutations/"
mkdir "$path_to_project"bed
mkdir "$path_to_project"annotated_somatic_mutations

#extracting somatic mutations 
grep ';SOMATIC;' "$path_to_vcf"variants.indel.vcf > "$path_to_vcf"variants.indel.filtered.vcf
grep ';SOMATIC;' "$path_to_vcf"variants.snp.vcf > "$path_to_vcf"variants.snp.filtered.vcf
#converting vcf to bed
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' "$path_to_vcf"variants.indel.filtered.vcf > "$path_to_bed"variants.indel.bed
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' "$path_to_vcf"variants.snp.filtered.vcf > "$path_to_bed"variants.snp.bed
#extracting all gene annotations for the somatic mutations
bedtools intersect -a "$path_to_gtf""$gtf" -b "$path_to_bed"variants.indel.bed > "$path_to_bed"variants.indel.intersect
bedtools intersect -a "$path_to_gtf""$gtf" -b "$path_to_bed"variants.snp.bed > "$path_to_bed"variants.snp.intersect
#extracting GTF lines corresponding to gene and gene names:
grep '\sgene\s' "$path_to_bed"variants.indel.intersect | awk '{print " " $1 " " $4 " " $5 " " $16}' > "$path_to_annotated_somatic_mutations"somatic_mutations.indel.txt
grep '\sgene\s' "$path_to_bed"variants.snp.intersect | awk '{print " " $1 " " $4 " " $5 " " $16}' > "$path_to_annotated_somatic_mutations"somatic_mutations.snp.txt

echo -e "\nDone! Annotated somatic mutations are saved to" "$path_to_annotated_somatic_mutations"
