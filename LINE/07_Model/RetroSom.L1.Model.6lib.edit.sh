#!/bin/bash -l
# NOTE the -l flag!
#
#SBATCH -J L1Model
#SBATCH --account=dflev
#SBATCH -p batch
#SBATCH --mem=10gb
#SBATCH --time=6:00:00

module load samtools
module load exonerate
module load bcftools
module load bedtools
module load perl-scg/5.14.4

Retro=$3
sub=$1\_Model
sub2=$1$8
echo $1  ### subject ID ###
echo $2  ### path ###
echo $Retro ### output folder ###
echo $4  ### reference genome hg19 G37 hg38... ###
echo $5  ### filter level:0=no filter, 1=prefilter, 2=allfilters ### 
masterpath=$7 
TE=LINE

mkdir $2/$sub
#rm -r $2/$sub/$Retro
mkdir $2/$sub/$Retro
mkdir $2/$sub/$Retro/$TE

### combine 6 invidual libraries ###
### combine SR support reads ###
awk '$5="lib1_"$5' OFS='\t' $2/$sub2\1/$Retro/$TE/$sub2\1.sr.pred.posreads > $2/$sub/$Retro/sr1.txt
awk '$5="lib2_"$5' OFS='\t' $2/$sub2\2/$Retro/$TE/$sub2\2.sr.pred.posreads > $2/$sub/$Retro/sr2.txt
awk '$5="lib3_"$5' OFS='\t' $2/$sub2\3/$Retro/$TE/$sub2\3.sr.pred.posreads > $2/$sub/$Retro/sr3.txt
awk '$5="lib4_"$5' OFS='\t' $2/$sub2\4/$Retro/$TE/$sub2\4.sr.pred.posreads > $2/$sub/$Retro/sr4.txt
awk '$5="lib5_"$5' OFS='\t' $2/$sub2\5/$Retro/$TE/$sub2\5.sr.pred.posreads > $2/$sub/$Retro/sr5.txt
awk '$5="lib6_"$5' OFS='\t' $2/$sub2\6/$Retro/$TE/$sub2\6.sr.pred.posreads > $2/$sub/$Retro/sr6.txt

cat $2/$sub/$Retro/sr1.txt \
    $2/$sub/$Retro/sr2.txt \
    $2/$sub/$Retro/sr3.txt \
    $2/$sub/$Retro/sr4.txt \
    $2/$sub/$Retro/sr5.txt \
    $2/$sub/$Retro/sr6.txt \
    > $2/$sub/$Retro/srall.txt 
cp $2/$sub/$Retro/srall.txt $2/$sub/$Retro/$sub.$TE.sr.tabe.discover
rm $2/$sub/$Retro/sr*.txt

### combine the PE reads ###
awk '$5="lib1_"$5' OFS='\t' $2/$sub2\1/$Retro/$TE/$sub2\1.pe.pred.posreads > $2/$sub/$Retro/pe1.txt
awk '$5="lib2_"$5' OFS='\t' $2/$sub2\2/$Retro/$TE/$sub2\2.pe.pred.posreads > $2/$sub/$Retro/pe2.txt
awk '$5="lib3_"$5' OFS='\t' $2/$sub2\3/$Retro/$TE/$sub2\3.pe.pred.posreads > $2/$sub/$Retro/pe3.txt
awk '$5="lib4_"$5' OFS='\t' $2/$sub2\4/$Retro/$TE/$sub2\4.pe.pred.posreads > $2/$sub/$Retro/pe4.txt
awk '$5="lib5_"$5' OFS='\t' $2/$sub2\5/$Retro/$TE/$sub2\5.pe.pred.posreads > $2/$sub/$Retro/pe5.txt
awk '$5="lib6_"$5' OFS='\t' $2/$sub2\6/$Retro/$TE/$sub2\6.pe.pred.posreads > $2/$sub/$Retro/pe6.txt

cat $2/$sub/$Retro/pe1.txt \
    $2/$sub/$Retro/pe2.txt \
    $2/$sub/$Retro/pe3.txt \
    $2/$sub/$Retro/pe4.txt \
    $2/$sub/$Retro/pe5.txt \
    $2/$sub/$Retro/pe6.txt \
    > $2/$sub/$Retro/$TE/$sub.$TE.novel.sites

rm $2/$sub/$Retro/pe1.txt
rm $2/$sub/$Retro/pe2.txt
rm $2/$sub/$Retro/pe3.txt
rm $2/$sub/$Retro/pe4.txt
rm $2/$sub/$Retro/pe5.txt
rm $2/$sub/$Retro/pe6.txt

awk '$5="lib1_"$5' OFS='\t' $2/$sub2\1/$Retro/$sub2\1.alignfilter.discover | grep L1 >  $2/$sub/$Retro/$sub.LINE.alignfilter.discover
awk '$5="lib2_"$5' OFS='\t' $2/$sub2\2/$Retro/$sub2\2.alignfilter.discover | grep L1 >> $2/$sub/$Retro/$sub.LINE.alignfilter.discover
awk '$5="lib3_"$5' OFS='\t' $2/$sub2\3/$Retro/$sub2\3.alignfilter.discover | grep L1 >> $2/$sub/$Retro/$sub.LINE.alignfilter.discover
awk '$5="lib4_"$5' OFS='\t' $2/$sub2\4/$Retro/$sub2\4.alignfilter.discover | grep L1 >> $2/$sub/$Retro/$sub.LINE.alignfilter.discover
awk '$5="lib5_"$5' OFS='\t' $2/$sub2\5/$Retro/$sub2\5.alignfilter.discover | grep L1 >> $2/$sub/$Retro/$sub.LINE.alignfilter.discover
awk '$5="lib6_"$5' OFS='\t' $2/$sub2\6/$Retro/$sub2\6.alignfilter.discover | grep L1 >> $2/$sub/$Retro/$sub.LINE.alignfilter.discover

### make calls ###
$masterpath/utls/25_filter_SR_dup.pl $1 $2/$sub $Retro $TE
grep L1 $2/$sub/$Retro/$sub.$TE.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u | \
   sort -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $2/$sub/$Retro/$TE/$sub.$TE.SR.calls

### PE calling ###
### L1 PE calls ###
$masterpath/utls/06_ana_depth.pl $sub $2/$sub $Retro $TE
$masterpath/utls/27_filter_dup.pl $1 $Retro $TE $2/$sub $5
#mv $2/$sub/$Retro/$TE/$sub.$TE.PE.calls $2/$sub/$Retro/$TE/$sub.$TE.PE.nodup.calls
$masterpath/utls/08_refine_depth.pl $sub $2/$sub $Retro $TE

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $2/$sub $Retro $TE $masterpath

windowBed \
   -w 100 -v \
   -a $2/$sub/$Retro/$TE/$sub.$TE.SR.PE.calls \
   -b $masterpath/refTE/position/$4.fa_LINE1_$6.bed \
   | windowBed -w 10 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.mask.bed \
   > $2/$sub/$Retro/$TE/$sub.$TE.noref.calls

### checking overlap ###
echo "checking overlapping reads"
$masterpath/utls/17_parse_overlap.pl $3 LINE $2 $1 $masterpath
### Read1/read2 from the same supporting read can only be used once ###
$masterpath/utls/16_filter_r1r2.pl $sub $2/$sub $Retro $TE

