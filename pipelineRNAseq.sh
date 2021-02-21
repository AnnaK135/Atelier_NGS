#!/bin/bash
# RNA-seq analysis 
# script 2: pipelineRNAseq.sh   
#--------------------------------------------------------------------------
# This script does the following:
# Trims fastq files to eliminate low quality bases at read extremities
# Maps the trimmed reads to the indexed reference genome -> .bam
# Performs BAM indexing (for visualisation in IGV)
# Counts reads falling to different genomic features 
# Annotates the counts table 

#Before running this script, please, run the script 1: preparationRNAseq.sh
#--------------------------------------------------------------------------

###For commenting out pieces of the code
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

path_to_project="/home/ubuntu/data/mydatalocal/TP1_RNAseq/"
path_to_raw_fastq="$path_to_project""raw_fastq/"
path_to_ref_genome="$path_to_project""ref_genome/"
path_to_index="$path_to_project""STARindex/"
path_to_bam="$path_to_project""aligned_reads/"
path_to_counts="$path_to_project""counts/"
annotation="annotation.gtf"
genome="chr18.fa"


###Sanity checks
echo "Checking for the necessary folders"
if [ ! -d "$path_to_project" ];then
  echo "Necessary folders do not exist! Please, run the script preparationRNAseq.sh first"
  exit 1
fi

if [ ! -d "$path_to_raw_fastq" ];then
  echo "Necessary folders do not exist! Please, run the script preparationRNAseq.sh first"
  exit 1
fi

if [ ! -d "$path_to_index"];then
  echo "STAR index does not exist. Please, run indexingRNAseq.sh"
  exit 1
fi
echo "Folders are OK"

###Trimming low quality bases 
echo "Trimming fastq files..."
ls "$path_to_raw_fastq" | sed -e 's/\.R[1,2].fastq$//' | uniq > "$path_to_project"fastq_prefixes
path_to_trimmed_fastq="$path_to_project""trimmed_fastq/"
mkdir "$path_to_trimmed_fastq"

for i in $(cat "$path_to_project"fastq_prefixes); do
  trimmomatic PE "$path_to_raw_fastq"${i}.R1.fastq "$path_to_raw_fastq"${i}.R2.fastq -baseout "$path_to_trimmed_fastq"${i}.trimmed.fastq LEADING:20 TRAILING:20 MINLEN:50
done

rm "$path_to_trimmed_fastq"*U.fastq #removing the files with trimmed reads
echo "Trimming done! Clean fastq files saved to" "$path_to_trimmed_fastq"

if [ ! -d "$path_to_trimmed_fastq" ];then
	echo "Trimming failed. Exiting"
	exit 1
fi

#### STAR Alignment
echo "Aligning reads to the genome..."
mkdir "$path_to_bam"

for i in $(cat "$path_to_project"fastq_prefixes);do

 STAR --runThreadN 4 --outFilterMultimapNmax 1\
  --genomeDir "$path_to_index" \
  --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$path_to_bam"${i} \
  --readFilesIn "$path_to_trimmed_fastq"${i}.trimmed_1P.fastq  "$path_to_trimmed_fastq"${i}.trimmed_2P.fastq

echo "Alignment for "${i}" completed! BAM file is saved to" "$path_to_bam"

### BAM indexing
echo "IndexiNG BAM..."
  samtools index "$path_to_bam"${i}Aligned.sortedByCoord.out.bam > "$path_to_bam"${i}Aligned.sortedByCoord.out.bam.bai
done
echo "BAM indexing completed! BAI files are saved to" "$path_to_bam"

### Counting reads
echo "Counting aligned reads..."
mkdir "$path_to_counts"
featureCounts -p -t exon -g gene_id -a "$path_to_ref_genome$annotation" -o "$path_to_counts"counts.txt "$path_to_bam"*.bam
echo "Done! Read counts table is saved in ""$path_to_counts""counts.txt" 

#Making read counts table reader-friendly
echo "Creating a tidy counts table with gene names..."
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
   "$path_to_ref_genome$annotation" | sort | uniq > "$path_to_project"encode-to-hugo.tab

sort "$path_to_counts"counts.txt > "$path_to_counts"temp1
sort "$path_to_project"encode-to-hugo.tab > "$path_to_counts"temp2
join "$path_to_counts"temp1 "$path_to_counts"temp2 |grep "chr18" > "$path_to_counts"temp3

awk '{print $13, $7, $8, $9, $10, $11, $12}' "$path_to_counts"temp3 > "$path_to_counts"counts_genes.txt
rm "$path_to_counts"temp*
echo "Done! Reader-friendly counts table is saved in ""$path_to_counts""counts_genes.txt"
ls "$path_to_bam"*.bam | sed -e 's/\.sampled.*.bam$//' > "$path_to_counts"samples_order.txt
echo "The info about columns order in counts_gene.txt is in samples_order.txt"


