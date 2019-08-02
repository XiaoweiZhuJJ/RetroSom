#!/bin/bash

usage="$(basename "$0") [-h] [-i t f c d e r s p m n g] -- Visualization tool for MEI supporting reads, e.g.,
./RetroVis.sh -i neg_NoModel -t LINE -f L1HS -c chr3 -d 153168189 -e 153168814 -r 61 -s 0 -p /home/xwzhu/transfer/BulkSeq/12004/ -m /home/xwzhu/masterpath -n test -g hg38

where:
    -h  show this help text
    -i  subject ID
    -t  TE class: LINE/ALU
    -f  TE family: L1HS/AluYa5
    -c  MEI chromosome
    -d  MEI coordinate 1
    -e  MEI coordinate 2
    -r  version control for RetroSom (default 1)
    -s  strandness of the MEI (+ strand: 1, - strand: 0)
    -p  datapath to the subject
    -m  masterpath (default ~/masterpath)
    -n  path to save the images
    -g  huamn reference genome (hg19/hg38)
"
masterpath='~/masterpath'
ver=1
while getopts ":hi:t:f:c:d:e:r:s:p:m:n:g:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    i) subject="$OPTARG"
       ;;
    t) TEclass="$OPTARG"
       ;;
    f) TEfamily="$OPTARG"
       ;;
    c) chr="$OPTARG"
       ;;
    d) cord1="$OPTARG"
       ;;
    e) cord2="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    p) datapath="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    n) pngpath="$OPTARG"
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

rm -r $masterpath/visual/temp
mkdir $masterpath/visual/temp
mkdir $masterpath/visual/$pngpath
echo -e $chr"\t"$cord1"\t"$cord2 > $masterpath/visual/temp/00_cord.bed
retro=retro_v$ver\_$strand
ref=$masterpath/refTE/sequence/$TEfamily.fa
callfile=$datapath/$subject/retro_v$ver\_$strand/$TEclass/$subject.$TEclass.SR.PE.calls
#callfile=$datapath/$subject/retro_v$ver\_$strand/$TEclass/$subject.$TEclass.noref.calls
PEfile=$datapath/$subject/retro_v$ver\_$strand/$TEclass/$subject.$TEclass.novel.sites
SRfile=$datapath/$subject/retro_v$ver\_$strand/$subject.$TEclass.sr.tabe.discover
#SRfile=$datapath/$subject/retro_v$ver\_$strand/$subject.sr.tabe.discover

module load perl
module load samtools
module load bedtools

### Step1: extract the MEI overlapping with the query insertion ###
awk '{printf $1"\t"$2"\t"$3"\t"; for (i=6; i<=NF; ++i) {printf "%s,", $i; if ($i == $NF) printf "\n"}}' \
  $callfile | \
  intersectBed -wo \
   -a stdin \
   -b $masterpath/visual/temp/00_cord.bed \
   > $masterpath/visual/temp/01_overlapping_insertions.txt

### Step2: extract the supporting reads ###
$masterpath/visual/MEI_support_reads.pl $PEfile $SRfile $masterpath $TEclass 

### Step3: Realign the reads to TEfamily consensus sequence ###
petemp=$masterpath/visual/temp/02_PEsupport.fa
if [ -s "$petemp" ]
then
$masterpath/exonerate/exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $petemp \
   $ref | grep "^INFO" \
   > $masterpath/visual/temp/04_PE.alignment
fi

srtemp=$masterpath/visual/temp/03_SRsupport.fa
if [ -s "$srtemp" ]
then
$masterpath/exonerate/exonerate --model affine:local \
   --bestn 1 --ryo "INFO: %qi %qal %pi %tS %ti %qab %qae %tab %tae %tas\n" \
   $srtemp \
   $ref | grep "^INFO" \
   > $masterpath/visual/temp/05_SR.alignment
fi

### step4: plotting the supporting reads ### 
imgfile=$subject\_retro$ver\_strand$strand\_$TEfamily\_$chr\_$cord1.svg
#echo $imgfile
$masterpath/visual/plot_MEI.SVG.pl $masterpath $imgfile $chr $TEfamily $pngpath $hg $strand 
