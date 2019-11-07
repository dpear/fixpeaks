#!/bin/bash
# run this from ~/Projects/fixpeaks$ ./scripts/p..sh

files=($(ls sums-peak-by-gene/*))

for f in ${files[@]}; do
	sample=$(basename $f | grep -Po 'q?[CLN][0-9]_OD[0-9]|WCE_OD[0-9]|K_OD[0-9]')
	#echo 'track' name=$sample 'useScore=1' > enriched/${sample}_enriched.bed
	#echo 'chrom chromStart chromEnd name score strand' >> enriched/${sample}_enriched.bed
	Rscript scripts/cluster3.R $f tmp
	cat tmp > sums3-enriched/${sample}_3_enriched.bed
done

rm tmp
