#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J L1PE
#SBATCH -o %x.%A.output
#SBATCH -e %x.%A.output
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --mem=20gb
#SBATCH --time=160:00:00

usage="$(basename "$0") [-h] [-o i m r g s] -- Predict the confidence of each supporting read to be real MEI using RF, NB and LR

where:
    -h  show this help text
    -o  output folder path
    -m  masterpath
    -r  RetroSom version control"

while getopts ":ho:m:r:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
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
echo "$outpath ran on $(hostname -s)"

module load r
Rscript $masterpath/LINE/02_PE_level1/01.RF.r $outpath $ver $masterpath
$masterpath/LINE/02_PE_level1/02_combine_all.sh $outpath $ver $masterpath
