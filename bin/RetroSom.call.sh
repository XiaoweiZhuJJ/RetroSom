module load samtools
module load exonerate
module load bcftools
module load bedtools

Retro=$3
strand=$6
echo $1  ### subject ID ###
echo $2  ### path ###
echo $Retro ### output folder ###
echo $4  ### reference genome hg19 G37 hg38... ###
echo $5  ### filter level:0=no filter, 1=prefilter, 2=allfilters ### 
echo $6  ### 1: + insert in forward strand; -1: - insert in reverse strand ###

masterpath="/home/xwzhu/levinson/RetroSom"
#echo before comment
#: <<'END'
date '+%m/%d/%y %H:%M:%S'
echo "TE calling phase ... begins"
### get the SEG scores ###
### add cigar to read names ###
awk -v str=$strand '{if (xor(str,($14>$15))) print ">"$5">"$4">"$7">"$12">"$13">"$14">"$15"\n"$16}' $2/$Retro/$1.sr.discover \
   > $2/$Retro/seg.sr.fasta

$masterpath/SEG/seg $2/$Retro/seg.sr.fasta 350 -h | \
   awk '{if (/>/) print $1"\t"$2}' | \
   sed s/complexity=// > $2/$Retro/seg.sr.scores

date '+%m/%d/%y %H:%M:%S'
echo "Modifying SR mode calls"
### filter SR mode ###
$masterpath/utls/05_retroseq_SR_filter_srlength.pl $1 $2 $Retro $5
$masterpath/utls/15_filter_SR_dup.pl $1 $2 $Retro

grep Alu $2/$Retro/$1.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u | \
   sort -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $2/$Retro/ALU/$1.ALU.SR.calls

grep L1 $2/$Retro/$1.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u | \
   sort -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $2/$Retro/LINE/$1.LINE.SR.calls

TE="HERVK"
grep $TE $2/$Retro/$1.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u | \
   sort -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $2/$Retro/$TE/$1.$TE.SR.calls

TE="HERVH"
grep $TE $2/$Retro/$1.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u | \
   sort -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $2/$Retro/$TE/$1.$TE.SR.calls

TE="SVA"
grep $TE $2/$Retro/$1.sr.dedup.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -u | \
   sort -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $2/$Retro/$TE/$1.$TE.SR.calls

### PE calling ###
### 1. the candidate reads ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... filter the candidate support"
$masterpath/utls/03_filter_PE_alignment.pl $1 $2 $Retro $5 $strand

$masterpath/SEG/seg $2/$Retro/seg.pe.fasta 350 -h | \
    awk '{if (/>/) print $1"\t"$2}' | \
    sed s/complexity=// > $2/$Retro/seg.pe.scores

$masterpath/utls/04_filter_PE_complexity.pl $1 $2 $Retro $5

date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... getting the candidate reads"
grep Alu $2/$Retro/$1.filter.discover | sort -k1 -n -k2 \
    > $2/$Retro/ALU/$1.ALU.novel.sites

grep L1 $2/$Retro/$1.filter.discover | sort -k1 -n -k2 \
    > $2/$Retro/LINE/$1.LINE.novel.sites

grep HERVK $2/$Retro/$1.filter.discover | sort -k1 -n -k2 \
    > $2/$Retro/HERVK/$1.HERVK.novel.sites

grep HERVH $2/$Retro/$1.filter.discover | sort -k1 -n -k2 \
    > $2/$Retro/HERVH/$1.HERVH.novel.sites

grep SVA $2/$Retro/$1.filter.discover | sort -k1 -n -k2 \
    > $2/$Retro/SVA/$1.SVA.novel.sites

### ALU PE calls ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Alu calls"
$masterpath/utls/06_ana_depth.pl $1 $2 $Retro ALU
$masterpath/utls/07_filter_dup.pl $1 $Retro ALU $2 $5
$masterpath/utls/08_refine_depth.pl $1 $2 $Retro ALU

### L1 PE calls ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... L1 calls"
$masterpath/utls/06_ana_depth.pl $1 $2 $Retro LINE
$masterpath/utls/07_filter_dup.pl $1 $Retro LINE $2 $5
$masterpath/utls/08_refine_depth.pl $1 $2 $Retro LINE

### HERV PE calls ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... HERV calls"
$masterpath/utls/06_ana_depth.pl $1 $2 $Retro HERVK
$masterpath/utls/07_filter_dup.pl $1 $Retro HERVK $2 $5
$masterpath/utls/08_refine_depth.pl $1 $2 $Retro HERVK

$masterpath/utls/06_ana_depth.pl $1 $2 $Retro HERVH
$masterpath/utls/07_filter_dup.pl $1 $Retro HERVH $2 $5
$masterpath/utls/08_refine_depth.pl $1 $2 $Retro HERVH

### SVA PE calls ###
$masterpath/utls/06_ana_depth.pl $1 $2 $Retro SVA
$masterpath/utls/07_filter_dup.pl $1 $Retro SVA $2 $5
$masterpath/utls/08_refine_depth.pl $1 $2 $Retro SVA

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $1 $2 $Retro ALU
$masterpath/utls/09_merge.SR.PE.support.sh $1 $2 $Retro LINE
$masterpath/utls/09_merge.SR.PE.support.sh $1 $2 $Retro HERVK
$masterpath/utls/09_merge.SR.PE.support.sh $1 $2 $Retro HERVH
$masterpath/utls/09_merge.SR.PE.support.sh $1 $2 $Retro SVA

#blurfl
#END
#echo after comment

### filter out the ref calls and the calls near reference masks ###
### filter out the insertions that are close to segdup, gaps(including telomere) and centromeres ##
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... ref TE filter"
windowBed \
   -w 100 -v \
   -a $2/$Retro/ALU/$1.ALU.SR.PE.calls \
   -b $masterpath/refTE/position/$4.fa_ALU\_$strand.bed \
   | windowBed -w 10 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.mask.bed \
   > $2/$Retro/ALU/$1.ALU.noref.calls

windowBed \
   -w 100 -v \
   -a $2/$Retro/LINE/$1.LINE.SR.PE.calls \
   -b $masterpath/refTE/position/$4.fa_LINE1\_$strand.bed \
   | windowBed -w 10 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.mask.bed \
   > $2/$Retro/LINE/$1.LINE.noref.calls

windowBed \
   -w 100 -v \
   -a $2/$Retro/HERVK/$1.HERVK.SR.PE.calls \
   -b $masterpath/refTE/position/$4.fa_HERV.bed \
   | windowBed -w 10 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.mask.bed \
   > $2/$Retro/HERVK/$1.HERVK.noref.calls

windowBed \
   -w 100 -v \
   -a $2/$Retro/HERVH/$1.HERVH.SR.PE.calls \
   -b $masterpath/refTE/position/$4.fa_HERV.bed \
   | windowBed -w 10 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.mask.bed \
   > $2/$Retro/HERVH/$1.HERVH.noref.calls

windowBed \
   -w 100 -v \
   -a $2/$Retro/SVA/$1.SVA.SR.PE.calls \
   -b $masterpath/refTE/position/$4.fa_SVA.bed | \
   windowBed \
   -w 100 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.fa_ALU.bed \
   | windowBed -w 10 -v \
   -a stdin \
   -b $masterpath/refTE/position/$4.mask.bed \
   > $2/$Retro/SVA/$1.SVA.noref.calls

### filter r1r2 ###
$masterpath/utls/16_filter_r1r2.pl $1 $2 $Retro ALU
$masterpath/utls/16_filter_r1r2.pl $1 $2 $Retro LINE
$masterpath/utls/16_filter_r1r2.pl $1 $2 $Retro HERVK
$masterpath/utls/16_filter_r1r2.pl $1 $2 $Retro HERVH
$masterpath/utls/16_filter_r1r2.pl $1 $2 $Retro SVA

### combine calls ###
cat $2/$Retro/ALU/$1.ALU.novel.calls \
    $2/$Retro/LINE/$1.LINE.novel.calls \
    $2/$Retro/HERVK/$1.HERVK.novel.calls \
    $2/$Retro/HERVH/$1.HERVH.novel.calls \
    $2/$Retro/SVA/$1.SVA.novel.calls \
    > $2/$Retro/MERGE/$1.novel.calls

cat $2/$Retro/ALU/$1.ALU.SR.PE.calls \
    $2/$Retro/LINE/$1.LINE.SR.PE.calls \
    $2/$Retro/HERVK/$1.HERVK.SR.PE.calls \
    $2/$Retro/HERVH/$1.HERVH.SR.PE.calls \
    $2/$Retro/SVA/$1.SVA.SR.PE.calls \
    > $2/$Retro/MERGE/$1.nodup.calls

### direction analysis ###
#$masterpath/utls/11_parse_direction.pl $1 ALU $Retro $2
#$masterpath/utls/11_parse_direction.pl $1 LINE $Retro $2

#blurfl
#END
#echo after comment
