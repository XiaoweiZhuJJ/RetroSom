#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i a] -- Remove PCR duplicates, secondary and supplementary alignment; Sort and indec the BAM file
Require BWA (v0.7.12), samtools and picard; Require 1 core

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -a  input bam file"

while getopts ":ho:i:a:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    a) bam="$OPTARG"
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

module load samtools
module load java
module load picard

### check the alignment statistics ###
samtools flagstat \
    $outpath/$sub/align/$sub.final.bam \
    > $outpath/$sub/QC/$sub.flagstat

### check the average insert size ###
picard CollectInsertSizeMetrics \
   I=$outpath/$sub/align/$sub.final.bam \
   O=$outpath/$sub/QC/$sub.insert.txt \
   H=$outpath/$sub/QC/$sub.histo.pdf \
   M=0.001

sed -n 8,8p $outpath/$sub/QC/$sub.insert.txt | awk '{print $1"\t"$2}' > $outpath/$sub/insert.$sub.txt

