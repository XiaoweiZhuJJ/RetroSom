#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J RemapL1
#SBATCH -o %x.%A.output
#SBATCH -e %x.%A.output
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --mem=20gb
#SBATCH --time=160:00:00

usage="$(basename "$0") [-h] [-o i m r g] -- Genotyping L1 supporting reads in several L1Hs unique loci. Requires one core, and 20gb memory

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath (default ~/masterpath)
    -r  RetroSom version control (default 1)
    -g  human reference genome (hg38 or hg19)"

ver=1
masterpath=~/masterpath
hg=hg38
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
echo "$sub ran on $(hostname -s)"

retro=retro_v$ver

Retro=fix_ACA
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/04_remap.NNN.sh $sub $outpath/$sub $Retro $hg ACA 5927 5928 5929 $masterpath

Retro=fix_ACG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/04_remap.NNN.sh $sub $outpath/$sub $Retro $hg ACG 5927 5928 5929 $masterpath

Retro=fix_TAG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/04_remap.NNN.sh $sub $outpath/$sub $Retro $hg TAG 6010 6011 6012 $masterpath

Retro=fix_5389
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/03_remap.NN.sh $sub $outpath/$sub $Retro $hg GC 5388 5389 $masterpath

Retro=fix_5533G
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/05_remap.NNNN.sh $sub $outpath/$sub $Retro $hg TGAG 5533 5534 5535 5536 $masterpath

Retro=fix_5533C
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/05_remap.NNNN.sh $sub $outpath/$sub $Retro $hg TGAC 5533 5534 5535 5536 $masterpath

Retro=fix_5710
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/retro_v$ver/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/retro_v$ver/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
$masterpath/LINE/00_remap/03_remap.NN.sh $sub $outpath/$sub $Retro $hg AT 5710 5711 $masterpath
