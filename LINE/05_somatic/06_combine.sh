module load bedtools
TE=LINE
sub=SOM$5\_$6
cut=$1
cont=$2
tissue=$3
control=$4
WIN=100
Retro=retro_v$5\_$6
workfolder=$7
#mkdir $workfolder/$sub
#mkdir $workfolder/$sub/$TE

awk -v cutoff=$cut '{if ($5>=cutoff) print}' \
   $workfolder/$tissue\_Model/$Retro/$TE/$tissue\_Model.$TE.novel.calls \
   > $workfolder/$tissue\_Model/$Retro/$TE/$tissue\_Model.$TE.$cut.calls

awk -v cutoff=$cont '{if ($5>=cutoff) print}' \
   $workfolder/$control\_NoModel/retro_v40_$6/$TE/$control\_NoModel.$TE.SR.PE.calls | \
   windowBed -v \
   -w $WIN \
   -a $workfolder/$tissue\_Model/$Retro/$TE/$tissue\_Model.$TE.$cut.calls \
   -b stdin \
   | sort -k1,1 -k2,3n \
   > $workfolder/$sub/$TE/$tissue.$TE.no$control.$cont\plus.learn.$cut.calls

