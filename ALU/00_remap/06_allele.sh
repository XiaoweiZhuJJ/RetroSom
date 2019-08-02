sub=$1
datapath=$2
Retro=$3
a1=$4
a2=$5

awk -v allele=$a1 '{if (($27 == "+") && ($28 == allele)) print}' $datapath/$Retro/$sub.Alu.discover \
   > $datapath/$Retro/$sub.discover
awk -v allale=$a2 '{if (($27 == "-") && ($28 == allele)) print}' $datapath/$Retro/$sub.Alu.discover \
   >> $datapath/$Retro/$sub.discover

awk -v allele=$a1 '{if (($27=="+") && ($28 == allele)) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   > $datapath/$Retro/$sub.sr.discover
awk -v allele=$a2 '{if (($27=="-") && ($28 == allele)) print}' $datapath/$Retro/$sub.sr.Alu.discover \
   >> $datapath/$Retro/$sub.sr.discover
