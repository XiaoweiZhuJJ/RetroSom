#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J Somatic
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
    -i  ID of the donor
    -m  masterpath
    -r  RetroSom version control"

while getopts ":ho:i:m:r:" opt; do
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

$masterpath/LINE/05_somatic/01_all_3tis_6lib.sh $ver 1 $outpath/$sub Astro Neuron Heart _
$masterpath/LINE/05_somatic/01_all_3tis_6lib.sh $ver 0 $outpath/$sub Astro Neuron Heart _
$masterpath/ALU/05_somatic/01_all_3tis_6lib.sh $ver 1 $outpath/$sub Astro Neuron Heart _
$masterpath/ALU/05_somatic/01_all_3tis_6lib.sh $ver 0 $outpath/$sub Astro Neuron Heart _

