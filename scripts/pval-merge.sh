#!/bin/bash

files=($(ls pval-peak-by-gene/*1_*))

for f1 in ${files[@]}; do
	sample=$(basename $f1 | grep -oP 'q?[CNL][0-9]_OD[0-9]')
	no_rep=$(echo $sample | sed -r 's/[0-9]//')	
	f2=$(echo $f1 | sed 's/1/2/')
	Rscript scripts/pval-merge.R $f1 $f2 ${no_rep}.bed
done
