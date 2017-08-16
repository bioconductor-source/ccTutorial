#!/bin/sh

# shell script to run Exonerate search

# 1st argument has to be the FASTA query file

FASTAFILE=$1
#echo "$FASTAFILE"
#OUTFILE=${1%.*}.ou  # substitute file ending
#echo "$OUTFILE"
#OUTFILE2=${1%.*}.psl

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M
do 
myfile="/nobackup/vol2/research/huber/genomes/mm9/chr$chr.fa"
echo $myfile
OUTFILE=chr${chr}.out
/nfs/research/huber/sw/bin/exonerate --percent 97 --showalignment no --showvulgar no --ryo ">\t%qi\t%ql\t%qS\t%ti\t%tab\t%tae\t%tS\t%ei\n" $1 $myfile > $OUTFILE
done
