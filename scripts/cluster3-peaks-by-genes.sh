#!/bin/bash
# run this from ~/Projects/fixpeaks$ ./scripts/p..sh

files=($(ls sums-peak-by-gene-14/*-b*))

for f in ${files[@]}; do
	sample=$(echo $f | grep -Po 'q?c?n?l?[CLN][0-9]-b|WCE[0-9]')
	Rscript scripts/cluster3.R $f tmp
	cat tmp > sums3-enriched-14/${sample}_3_enriched.bed
	echo $sample
done

rm tmp
