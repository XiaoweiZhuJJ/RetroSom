#!/bin/bash -l
# NOTE the -l flag!
#
#SBATCH -J RetroDisover
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=200:00:00

module load samtools
export PATH=$PATH:/srv/gsfs0/projects/levinson/xwzhu/exonerate/
module load bcftools
module load bedtools
module load perl
sub=$1
outpath=$2
Retro=$3
masterpath=$4
echo "$sub ran on $(hostname -s)"

date '+%m/%d/%y %H:%M:%S'
echo "Discovery phase ... begins"
$masterpath/bin/retroseq.prun.pl -discover \
  -align \
  -srmode \
  -minclip 20 \
  -len 26 \
  -srcands $outpath/$sub/$Retro/$sub.sr.discover \
  -bam $outpath/$sub/align/$sub.final.bam \
  -eref $masterpath/refTE/TE_ALHS.bed \
  -output $outpath/$sub/$Retro/$sub.discover

### seperate the supporting reads from either strand ###
$masterpath/utls/21_par_direction.pl $sub $outpath $Retro
date '+%m/%d/%y %H:%M:%S'
echo "Discovery phase ... finish"

