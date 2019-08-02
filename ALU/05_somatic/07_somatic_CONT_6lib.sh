#!/bin/bash -l
# NOTE the -l flag!
#
#SBATCH -J SOM
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --mem=12gb
#SBATCH --time=6:00:00

module load bedtools
TE=ALU
sub=SOM$5\_$6
cut=$1
cont=$2
tissue=$3
control=$4$8
WIN=600
Retro=retro_v$5\_$6
workfolder=$7

awk '{if (($4 ~/Alu/) && ($1~/chr/)) print}' \
   $workfolder/$control\1/retro_$5/$control\1.$6.discover | \
   windowBed -w $WIN \
   -a stdin \
   -b $workfolder/$sub/$TE/$tissue.$TE.no$4.$cont\plus.learn.$cut.calls \
   | awk '{if (($4 ~ /Alu/) && ($8 > 95) ) print cont1"\t"$0}' \
   > $workfolder/$sub/ALU/$tissue.ALU.no$4.$cont\plus.learn.$cut.old.discover

awk '{if (($4 ~/Alu/) && ($1~/chr/)) print}' \
   $workfolder/$control\2/retro_$5/$control\2.$6.discover | \
   windowBed -w $WIN \
   -a stdin \
   -b $workfolder/$sub/$TE/$tissue.$TE.no$4.$cont\plus.learn.$cut.calls \
   | awk '{if (($4 ~ /Alu/) && ($8 > 95) ) print cont2"\t"$0}' \
   >> $workfolder/$sub/ALU/$tissue.ALU.no$4.$cont\plus.learn.$cut.old.discover

awk '{if (($4 ~/Alu/) && ($1~/chr/)) print}' \
   $workfolder/$control\3/retro_$5/$control\3.$6.discover | \
   windowBed -w $WIN \
   -a stdin \
   -b $workfolder/$sub/$TE/$tissue.$TE.no$4.$cont\plus.learn.$cut.calls \
   | awk '{if (($4 ~ /Alu/) && ($8 > 95) ) print cont3"\t"$0}' \
   >> $workfolder/$sub/ALU/$tissue.ALU.no$4.$cont\plus.learn.$cut.old.discover

awk '{if (($4 ~/Alu/) && ($1~/chr/)) print}' \
   $workfolder/$control\4/retro_$5/$control\4.$6.discover | \
   windowBed -w $WIN \
   -a stdin \
   -b $workfolder/$sub/$TE/$tissue.$TE.no$4.$cont\plus.learn.$cut.calls \
   | awk '{if (($4 ~ /Alu/) && ($8 > 95) ) print cont4"\t"$0}' \
   >> $workfolder/$sub/ALU/$tissue.ALU.no$4.$cont\plus.learn.$cut.old.discover

awk '{if (($4 ~/Alu/) && ($1~/chr/)) print}' \
   $workfolder/$control\5/retro_$5/$control\5.$6.discover | \
   windowBed -w $WIN \
   -a stdin \
   -b $workfolder/$sub/$TE/$tissue.$TE.no$4.$cont\plus.learn.$cut.calls \
   | awk '{if (($4 ~ /Alu/) && ($8 > 95) ) print cont5"\t"$0}' \
   >> $workfolder/$sub/ALU/$tissue.ALU.no$4.$cont\plus.learn.$cut.old.discover

awk '{if (($4 ~/Alu/) && ($1~/chr/)) print}' \
   $workfolder/$control\6/retro_$5/$control\6.$6.discover | \
   windowBed -w $WIN \
   -a stdin \
   -b $workfolder/$sub/$TE/$tissue.$TE.no$4.$cont\plus.learn.$cut.calls \
   | awk '{if (($4 ~ /Alu/) && ($8 > 95) ) print cont6"\t"$0}' \
   >> $workfolder/$sub/ALU/$tissue.ALU.no$4.$cont\plus.learn.$cut.old.discover

