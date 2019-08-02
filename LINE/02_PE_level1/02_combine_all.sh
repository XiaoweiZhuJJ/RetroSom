outpath=$1
ver=$2
masterpath=$3
while read line; do
   $masterpath/LINE/02_PE_level1/03_combine_PE.pl $line $ver $outpath
done < $outpath/list.txt

