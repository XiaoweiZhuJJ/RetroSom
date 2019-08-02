#!/bin/bash -l

usage="$(basename "$0") [-h] [-1 -2 -o -i] -- program to trim illumina adapters and low sequencing quality tails

where:
    -h  show this help text
    -1  input sequencing file for read1
    -2  input sequencing file for read2
    -o  output folder path
    -i  subject ID"

while getopts ":h1:2:o:i:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    1) in_1="$OPTARG"
       ;;
    2) in_2="$OPTARG"
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
### require java, python, cutadapt, and Fastqc ###
module load java 
module load python
module load fastqc
module load cutadapt
module load samtools

### trim sequencing adapters, low sequencing quality tails (<20), and only keep the reads longer than 30 bp ###
cutadapt \
   -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
   --trim-n -q 20 -m 30 \
   -o $outpath/$sub/reads/$sub.trimmed.R1.fastq \
   -p $outpath/$sub/reads/$sub.trimmed.R2.fastq \
   --info-file $outpath/$sub/QC/trimming.info \
   $in_1 \
   $in_2 

### checking the reads quality with FastQC ###
fastqc -f fastq -t 4 \
   -o $outpath/$sub/QC \
   $outpath/$sub/reads/$sub.trimmed.R1.fastq

fastqc -f fastq -t 4 \
   -o $outpath/$sub/QC \
   $outpath/$sub/reads/$sub.trimmed.R2.fastq

## Alignment and post processing ###
### using BWA version 0.7.12 ###
### using human reference genome GRCh38DH ###
bwa mem -t 6 -B 4 -O 6 -E 1 -M \
   -R $Itag \
   $refG/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   $outpath/$sub/reads/$sub.trimmed.R1.fastq \
   $outpath/$sub/reads/$sub.trimmed.R2.fastq \
   | /scg/apps/software/bwakit/0.7.15/bin/k8 /scg/apps/software/bwa/0.7.12/bwakit/bwa-postalt.js -p \
   $refG/GRCh38_full_analysis_set_plus_decoy_hla-extra.fa \
   $refG/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt \
   | samtools view -bSh \
   > $outpath/$sub/align/$sub.bam

### sort the aligned file by coordinates ###
samtools sort -@ 6 -T $sub \
   $outpath/$sub/align/$sub.bam \
   -o $outpath/$sub/align/$sub.sort.bam
