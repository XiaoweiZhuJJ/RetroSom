#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J Preprocess
#SBATCH -o %x.%A.output
#SBATCH -e %x.%A.output
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --mem=20gb
#SBATCH --time=20:00:00

usage="$(basename "$0") [-h] [-o i m r g] -- Preprocessing before calling the putative MEIs.
Requires 1 core, and 20gb memory; Requires samtools, bedtools, bcftools and exonerate

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath (default ~/masterpath)
    -r  RetroSom version control (default 1)
    -g  reference genome version (hg38, hg19 etc., default: hg38)"

ver=1
masterpath=~/masterpath
while getopts ":ho:i:m:r:g:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    g) hg="$OPTARG"
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
retro=retro_v$ver
workfolder=$outpath/$sub
echo "$sub ran on $(hostname -s)"

module load samtools
module load bcftools
module load bedtools

date '+%m/%d/%y %H:%M:%S'
echo "TE calling phase ... begins"

### Group1 => 5, convert part of PE supporting reads to SR ###
$masterpath/utls/00_GroupItoV.pl $masterpath $sub $workfolder/$retro
### select a window around supporting reads ###
$masterpath/utls/01_correct.bed.pl $masterpath $sub 300 $workfolder $retro $hg
### check the sequencing depth at the selected windows, require bedtools ### 
$masterpath/utls/02_depth_zero.sh $sub $workfolder $retro
date '+%m/%d/%y %H:%M:%S'
echo "Depth analysis ... finishes"
