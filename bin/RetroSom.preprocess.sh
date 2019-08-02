module load samtools
module load exonerate
module load bcftools
module load bedtools

sub=$1          ### subject ID ###
workfolder=$2   ### path ###
Retro=$3        ### output folder ###
hg=$4           ### reference genome hg19 G37 hg38... ###

masterpath="/home/xwzhu/levinson/RetroSom"
date '+%m/%d/%y %H:%M:%S'
echo "TE calling phase ... begins"

### Group1 => 5, convert part of PE supporting reads to SR ###
$masterpath/utls/00_GroupItoV.pl $sub $workfolder/$Retro

### calculate the depth around support reads ###
$masterpath/utls/01_correct.bed.pl $sub 300 $workfolder $Retro $hg
$masterpath/utls/02_depth_zero.sh $sub $workfolder $Retro
date '+%m/%d/%y %H:%M:%S'
echo "Depth analysis ... finishes"

