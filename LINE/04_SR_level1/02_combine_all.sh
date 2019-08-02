outpath=$1
ver=$2
masterpath=$3
while read line; do
    $masterpath/LINE/04_SR_level1/03_combine_SR.pl $line $ver $outpath
done < $outpath/list.txt

