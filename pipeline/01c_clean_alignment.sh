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

### remove PCR duplicates, secondary alignment and supplementary alignment ###
picard MarkDuplicates \
   I=$bam \
   O=/dev/stdout \
   ASSUME_SORTED=true \
   TMP_DIR=$outpath/$sub/align/ \
   MAX_RECORDS_IN_RAM=2000000 \
   REMOVE_DUPLICATES=true \
   VALIDATION_STRINGENCY=SILENT \
   METRICS_FILE=$outpath/$sub/QC/$sub.dedupp.metrics | \
   samtools view -bh -F 0x0100 -F 0x400 -F 0x800 \
   > $outpath/$sub/align/$sub.final.bam

### index the BAM file ###
samtools index \
   $outpath/$sub/align/$sub.final.bam

