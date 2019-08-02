TE=LINE
sub=SOM$1\_$2
Retro=retro_v$1\_$2
workfolder=$3
astro="$4$9"
neuron="$5$9"
heart="$6$9"
hg=$7
masterpath=$8
echo "$astro $neuron $heart"
mkdir $workfolder/$sub
mkdir $workfolder/$sub/$TE

### extract supporting reads that passed the RF modeling ###
$masterpath/utls/12_postive_support.pl $astro\1 $workfolder/$astro\1 $Retro
$masterpath/utls/12_postive_support.pl $astro\2 $workfolder/$astro\2 $Retro
$masterpath/utls/12_postive_support.pl $astro\3 $workfolder/$astro\3 $Retro
$masterpath/utls/12_postive_support.pl $astro\4 $workfolder/$astro\4 $Retro
$masterpath/utls/12_postive_support.pl $astro\5 $workfolder/$astro\5 $Retro
$masterpath/utls/12_postive_support.pl $astro\6 $workfolder/$astro\6 $Retro

$masterpath/utls/12_postive_support.pl $neuron\1 $workfolder/$neuron\1 $Retro
$masterpath/utls/12_postive_support.pl $neuron\2 $workfolder/$neuron\2 $Retro
$masterpath/utls/12_postive_support.pl $neuron\3 $workfolder/$neuron\3 $Retro
$masterpath/utls/12_postive_support.pl $neuron\4 $workfolder/$neuron\4 $Retro
$masterpath/utls/12_postive_support.pl $neuron\5 $workfolder/$neuron\5 $Retro
$masterpath/utls/12_postive_support.pl $neuron\6 $workfolder/$neuron\6 $Retro

$masterpath/utls/12_postive_support.pl $heart\1 $workfolder/$heart\1 $Retro
$masterpath/utls/12_postive_support.pl $heart\2 $workfolder/$heart\2 $Retro
$masterpath/utls/12_postive_support.pl $heart\3 $workfolder/$heart\3 $Retro
$masterpath/utls/12_postive_support.pl $heart\4 $workfolder/$heart\4 $Retro
$masterpath/utls/12_postive_support.pl $heart\5 $workfolder/$heart\5 $Retro
$masterpath/utls/12_postive_support.pl $heart\6 $workfolder/$heart\6 $Retro
echo "calling combined TE"

$masterpath/LINE/07_Model/RetroSom.L1.Model.6lib.edit.sh $4 $workfolder $Retro $hg 1 $2 $masterpath $9
$masterpath/LINE/07_Model/RetroSom.L1.Model.6lib.edit.sh $5 $workfolder $Retro $hg 1 $2 $masterpath $9
$masterpath/LINE/07_Model/RetroSom.L1.Model.6lib.edit.sh $6 $workfolder $Retro $hg 1 $2 $masterpath $9
$masterpath/LINE/07_Model/RetroSom.L1.NoModel.6lib.sh $4 $workfolder $Retro $hg 1 $2 $masterpath $9
$masterpath/LINE/07_Model/RetroSom.L1.NoModel.6lib.sh $5 $workfolder $Retro $hg 1 $2 $masterpath $9
$masterpath/LINE/07_Model/RetroSom.L1.NoModel.6lib.sh $6 $workfolder $Retro $hg 1 $2 $masterpath $9

$masterpath/LINE/05_somatic/05_submit_somatic_6lib.sh $6 $4 $1 $2 $3 $masterpath $9
$masterpath/LINE/05_somatic/05_submit_somatic_6lib.sh $6 $5 $1 $2 $3 $masterpath $9
$masterpath/LINE/05_somatic/05_submit_somatic_6lib.sh $4 $6 $1 $2 $3 $masterpath $9
$masterpath/LINE/05_somatic/05_submit_somatic_6lib.sh $5 $6 $1 $2 $3 $masterpath $9
