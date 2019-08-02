#! /bin/sh
#
#$ -N ALUYB
#
#$ -cwd
#$ -l h_rt=6:00:00
#
#$ -j y
#
#$ -S /bin/bash
#

module load blat

sub=$1
datapath=$2
Retro=$3

#echo before comment
#: <<'END'
mv $datapath/$Retro/$sub.discover \
   $datapath/$Retro/$sub.alleles.discover
mv $datapath/$Retro/$sub.sr.discover \
   $datapath/$Retro/$sub.sr.alleles.discover

### ALU ###
awk '{if ($11<245 && $12>=260 && ($4 ~ /AluYb/)) print}' \
   $datapath/$Retro/$sub.alleles.discover \
   | awk '{print ">"$0"\tpe\n"$13}' | sed 's/\t/\_/g' \
   > $datapath/$Retro/$sub.ALU.fa

awk '{if ($11>=260 && $12<245 && ($4 ~ /AluYb/)) print}' \
   $datapath/$Retro/$sub.alleles.discover \
   | awk '{print ">"$0"\tpe\n"$13}' | sed 's/\t/\_/g' \
   >> $datapath/$Retro/$sub.ALU.fa

awk '{if ($14<245 && $15>=260 && ($4 ~ /AluYb/)) print}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   | awk '{print ">"$0"\tsr\n"$17}' | sed 's/\t/\_/g' \
   >> $datapath/$Retro/$sub.ALU.fa

awk '{if ($14>=260 && $15<245 && ($4 ~ /AluYb/)) print}' \
   $datapath/$Retro/$sub.sr.alleles.discover \
   | awk '{print ">"$0"\tsr\n"$17}' | sed 's/\t/\_/g' \
   >> $datapath/$Retro/$sub.ALU.fa

blat -out=blast9 \
   /home/xwzhu/levinson/retro_filter/retro_camal/refTE/sequence/AluYb.fa \
   $datapath/$Retro/$sub.ALU.fa \
   $datapath/$Retro/$sub.ALU.alignout

#/home/xwzhu/levinson/retro_camal2/Alu_remap/03_realign.ALU.pl $sub $Retro $datapath

#grep L1 $datapath/$Retro/$sub.alleles.discover \
#   >> $datapath/$Retro/$sub.discover
#grep L1 $datapath/$Retro/$sub.sr.alleles.discover \
#   >> $datapath/$Retro/$sub.sr.discover

