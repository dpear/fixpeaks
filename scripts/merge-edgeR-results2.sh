#!/bin/bash

#Usage: ./merge-edgeR-results2.sh <tf1 file> <tf2 file> <tf3 file> <output file name>
# where the tf files are the edgeR outputs from run-edgeR.R wrapper script

t1=$1
t2=$2
t3=$3
out=$4

awk '{print $0,FILENAME}' $t1 | sed 's/_OD1-edgeR.tsv//' > $out
awk '{print $0,FILENAME}' $t2 | sed 's/_OD1-edgeR.tsv//' >> $out 
awk '{print $0,FILENAME}' $t3 | sed 's/_OD1-edgeR.tsv//' >> $out
cat $out | grep -v 'gene' > tmp
echo $(head -1 $t1) tf > $out
cat $tmp >> $out
rm tmp


