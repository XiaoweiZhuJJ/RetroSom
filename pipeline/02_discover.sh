#!/bin/bash

usage="$(basename "$0") [-h] [-o i g] -- Discover putative MEI supporting reads, based on Retroseq; Requires one core, and 4-10G memory per core. 

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath (default ~/masterpath)
    -r  version control for RetroSom (default 1)

Please note that post-processing is slow, and it may take ~100hours to align sequencing reads from 30X W
GS data"
masterpath='~/masterpath'
ver=1
while getopts ":ho:i:m:r:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    r) ver="$OPTARG"
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

### make new folders if this is a new analysis ###
### delete $outpath/$sub/$retro if it already existed ###
retro=retro_v$ver
mkdir $outpath/$sub/script
mkdir $outpath/$sub/$retro
mkdir $outpath/$sub/$retro/ALU
mkdir $outpath/$sub/$retro/LINE
mkdir $outpath/$sub/$retro/MERGE

### copy script to folder ###
cp $masterpath/bin/RetroSom.discover.sh $outpath/$sub/script/
cd $outpath/$sub/script

### run the script to discover candidate supporting reads ###
sbatch ./RetroSom.discover.sh $sub $outpath $retro $masterpath
