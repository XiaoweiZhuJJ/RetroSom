#! /bin/bash

usage="$(basename "$0") [-h] [-o i m r g t a b c] -- Discovering somatic MEI insertions supported with >=2 reads.

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath (default ~/masterpath)
    -r  version control for RetroSom (default 1)
    -g  reference genome (supporting hg38, hg19 and b37)
    -t  type of input (1=raw_sequencing_reads;2=BAM_to_be_realigned; 3=Aligned_BAM)
    -a  seqeuncing read1 (required if input are raw sequencing reads, input==1)
    -b  sequencing read2 (required if input are raw sequencing reads, input==1)
    -c  input BAM file (input==2 or 3)

### RetroSom dependency ###
1. Perl (v5.28.1)
Packages: 'GD', 'GD::Arrow', 'GD::SVG', 'Parallel::ForkManager'
 
2. R (v3.5.0)
Packages: 'randomForest', 'glmnet',  'e1071', 'PRROC'
 
3. Genomics tools:
bedtools (https://github.com/arq5x/bedtools2)
picard (https://github.com/broadinstitute/picard/releases/tag/2.20.4)
samtools (http://samtools.sourceforge.net/)
exonerate (https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
SEG (http://www.biology.wustl.edu/gcg/seg.html)
 
4. Pipeline management
SLURM
"

while getopts ":ho:i:m:r:t:g:a:b:c:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    t) datatype="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    a) read1="$OPTARG"
       ;;
    b) read2="$OPTARG"
       ;;
    c) bam="$OPTARG"
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

######################################
### Step0: Creating output folders ###
######################################
mkdir $outpath/$sub
mkdir $outpath/$sub/reads
mkdir $outpath/$sub/align
mkdir $outpath/$sub/QC
mkdir $outpath/$sub/script
mkdir $outpath/$sub/retro_v$ver
mkdir $outpath/$sub/retro_v$ver/ALU
mkdir $outpath/$sub/retro_v$ver/LINE
cd $outpath/$sub/script

#################################
### Step1: Sequence Alignment ###
################################# 
# Input Type 1: raw sequencing reads ### 
# Input Type 2: aligned reads to be realigned ### 
# Input Type 3: aligned reads (no realignment) ### 
slurm_sc="-o %x.%A.output -e %x.%A.output -A dflev -p batch --mem=20gb --time=60:00:00"
slurm_mc="-o %x.%A.output -e %x.%A.output -A dflev -p batch --time=200:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=10gb"
if [ "$datatype" == 1 ]
then
    jid1a=$(sbatch -J AlignReads $slurm_sc $masterpath/pipeline/01a_align_reads.sh -1 $read1 -2 $read2 -o $outpath -i $sub |awk '{print $4}')
    jid1c=$(sbatch -J CleanBam $slurm_sc --dependency=afterok:$jid1a $masterpath/pipeline/01c_clean_alignment.sh -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub | awk '{print $4}')
    jid1d=$(sbatch -J QCBam $slurm_sc --dependency=afterok:$jid1c $masterpath/pipeline/01d_qc_alignment.sh -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub | awk '{print $4}')
elif [ "$datatype" == 2 ]
then
    jid1b=$(sbatch -J RealignBAM $slurm_sc $masterpath/pipeline/01b_realign_BAM.sh -a $bam -o $outpath -i $sub | awk '{print $4}')
    jid1c=$(sbatch -J CleanBAM $slurm_sc --dependency=afterok:$jid1b $masterpath/pipeline/01c_clean_alignment.sh -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub | awk '{print $4}')
    jid1d=$(sbatch -J QCBams $slurm_sc --dependency=afterok:$jid1c $masterpath/pipeline/01d_qc_alignment.sh -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub | awk '{print $4}')
elif [ "$datatype" == 3 ]
then
    jid1c=$(sbatch -J CleanBam $slurm_sc $masterpath/pipeline/01c_clean_alignment.sh -a $bam -o $outpath -i $sub | awk '{print $4}')
    jid1d=$(sbatch -J QCBam $slurm_sc --dependency=afterok:$jid1c $masterpath/pipeline/01d_qc_alignment.sh -a $bam -o $outpath -i $sub | awk '{print $4}')
elif [ "$datatype" == 4 ]
then
    jid1c=$(sbatch -J CleanBam $slurm_sc $masterpath/pipeline/01e_linkBAM.sh -a $bam -o $outpath -i $sub | awk '{print $4}')
    jid1d=$(sbatch -J QCBam $slurm_sc --dependency=afterok:$jid1c $masterpath/pipeline/01d_qc_alignment.sh -a $bam -o $outpath -i $sub | awk '{print $4}')
fi

##################################################
### Step2: Discover candidate supporting reads ###
##################################################
retro=retro_v$ver
jid2=$(sbatch -J RetroDiscover $slurm_mc --dependency=afterok:$jid1c $masterpath/bin/RetroSom.discover.sh $sub $outpath $retro $masterpath | awk '{print $4}')

##############################################
### Step3: Putative MEIs without filtering ###
##############################################
jid3a=$(sbatch -J PreProc $slurm_sc --dependency=afterok:$jid2 $masterpath/pipeline/03a_preprocess.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg | awk '{print $4}')
# insertions in -strand #
jid3b=$(sbatch -J PutMEI $slurm_sc --dependency=afterok:$jid2:$jid3a $masterpath/pipeline/03b_call_putative_MEI.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 -f 1 | awk '{print $4}')
# insertions in +strand #
jid3c=$(sbatch -J PutMEI $slurm_sc --dependency=afterok:$jid2:$jid3a $masterpath/pipeline/03b_call_putative_MEI.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 -f 1 | awk '{print $4}')

######################################################
### Step4: Remapping L1HS or AluY specific Alleles ###
######################################################
# Realign L1 supporting reads #
jid4a=$(sbatch -J L1MAP $slurm_sc --dependency=afterok:$jid2 $masterpath/pipeline/04a_remap_L1.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg | awk '{print $4}')
# Realign Alu supporting reads #
jid4b=$(sbatch -J AluMAP $slurm_sc --dependency=afterok:$jid2 $masterpath/pipeline/04b_remap_Alu.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg | awk '{print $4}')

####################################################
### Step5: Matricies for L1/Alu supporting reads ###
####################################################
# -strand #
jid5a0=$(sbatch -J L1PEMat $slurm_sc --dependency=afterok:$jid3b:$jid4a $masterpath/pipeline/05a_matrix_gen_L1PE.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 | awk '{print $4}')
jid5b0=$(sbatch -J L1SRMat $slurm_sc --dependency=afterok:$jid3b:$jid4a $masterpath/pipeline/05b_matrix_gen_L1SR.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 | awk '{print $4}')
jid5c0=$(sbatch -J AluPEMat $slurm_sc --dependency=afterok:$jid3b:$jid4b $masterpath/pipeline/05c_matrix_gen_AluPE.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 | awk '{print $4}')
jid5d0=$(sbatch -J AluSRMat $slurm_sc --dependency=afterok:$jid3b:$jid4b $masterpath/pipeline/05d_matrix_gen_AluSR.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 | awk '{print $4}')
# +strand #
jid5a1=$(sbatch -J L1PEMat $slurm_sc --dependency=afterok:$jid3c:$jid4a $masterpath/pipeline/05a_matrix_gen_L1PE.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 | awk '{print $4}')
jid5b1=$(sbatch -J L1SRMat $slurm_sc --dependency=afterok:$jid3c:$jid4a $masterpath/pipeline/05b_matrix_gen_L1SR.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 | awk '{print $4}')
jid5c1=$(sbatch -J AluPEMat $slurm_sc --dependency=afterok:$jid3c:$jid4b $masterpath/pipeline/05c_matrix_gen_AluPE.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 | awk '{print $4}')
jid5d1=$(sbatch -J AluSRMat $slurm_sc --dependency=afterok:$jid3c:$jid4b $masterpath/pipeline/05d_matrix_gen_AluSR.sh -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 | awk '{print $4}')

########################################################
### Step6: Classification with a random forest model ###
########################################################
jid6a=$(sbatch -J L1PEl1 $slurm_sc --dependency=afterok:$jid5a0:$jid5a1 $masterpath/pipeline/06a_Level1_L1PE.sh -o $outpath -m $masterpath -r $ver | awk '{print $4}')
jid6b=$(sbatch -J L1SRl1  $slurm_sc --dependency=afterok:$jid5b0:$jid5b1 $masterpath/pipeline/06b_Level1_L1SR.sh -o $outpath -m $masterpath -r $ver | awk '{print $4}')
jid6c=$(sbatch -J AluPEl1 $slurm_sc --dependency=afterok:$jid5c0:$jid5c1 $masterpath/pipeline/06c_Level1_AluPE.sh -o $outpath -m $masterpath -r $ver | awk '{print $4}')
jid6d=$(sbatch -J AluSRl1 $slurm_sc --dependency=afterok:$jid5d0:$jid5d1 $masterpath/pipeline/06d_Level1_AluSR.sh -o $outpath -m $masterpath -r $ver | awk '{print $4}')

##################################################
### Step7. Calling Putative somatic insertions ###
##################################################
jid7=$(sbatch -J Somatic $slurm_sc --dependency=afterok:$jid6a:$jid6b:$jid6c:$jid6d $masterpath/pipeline/07_somatic.sh -o $outpath -i $sub -m $masterpath -r $ver | awk '{print $4}')
