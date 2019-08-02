#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J RemapAlu
#SBATCH -o %x.%A.output
#SBATCH -e %x.%A.output
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --mem=20gb
#SBATCH --time=160:00:00

usage="$(basename "$0") [-h] [-o i m r g] -- Genotyping Alu supporting reads in several AluY unique loci. Requires one core, and 20gb memory

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
Retro=fix_GGCG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/01_GGCG.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_GAGC
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/02_GAGC.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_CC
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/03_CC.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_CG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/03_CG.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_TG
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/03_TG.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_AluYb
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/04_AluYb.sh $sub $outpath/$sub $Retro $masterpath

Retro=fix_AT
rm -r $outpath/$sub/$Retro
mkdir $outpath/$sub/$Retro
ln -s $outpath/$sub/$retro/$sub.discover $outpath/$sub/$Retro/$sub.discover
ln -s $outpath/$sub/$retro/$sub.sr.discover $outpath/$sub/$Retro/$sub.sr.discover
 $masterpath/ALU/00_remap/05_AT.sh $sub $outpath/$sub $Retro $masterpath

