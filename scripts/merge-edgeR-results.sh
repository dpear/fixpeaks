#!/bin/bash

# Usage: ./merge-edgeR-results.sh <tsv1> <tsv2> <tsv3> <output file name>

t1=$1
t2=$2
t3=$3
out=$4

echo $(cat $t1 | head -1 | sed -r "s/\s/\s $t1/g" | sed "s/edgeR.tsv//g") ${t1}_tmp > ${t1}_tmp
echo $(cat $t2 | head -1 | sed -r "s/\s/\s $t2/g" | sed "s/edgeR.tsv//g") ${t2}_tmp > ${t2}_tmp
echo $(cat $t3 | head -1 | sed -r "s/\s/\s $t3/g" | sed "s/edgeR.tsv//g") ${t3}_tmp > ${t3}_tmp
sort $t1 >> ${t1}_tmp; awk '{print $0,FILENAME}' ${t1}_tmp > t; mv t ${t1}_tmp
sort $t2 >> ${t2}_tmp; awk '{print $0,FILENAME}' ${t2}_tmp > t; mv t ${t2}_tmp
sort $t3 >> ${t3}_tmp; awk '{print $0,FILENAME}' ${t3}_tmp > t; mv t ${t3}_tmp

join ${t1}_tmp ${t2}_tmp -j 1 --header > j-tmp
join j-tmp ${t3}_tmp -j 1 --header > $out

# get rid of the stray unlabelled header
grep -vP ' FDR' $out > tmp; mv tmp $out

rm ${t1}_tmp
rm ${t2}_tmp
rm ${t3}_tmp
