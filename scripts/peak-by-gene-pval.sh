#!/bin/bash
macs=($(ls MACS_K_control/*/*.xls))

for m in ${macs[@]}; do
	sample=($(echo $m | grep -Po 'q?[CNL][0-9]_OD[0-9]'))
	sample=${sample[1]}
	reg=peak-by-gene-outputs/${sample}_final.bed
	Rscript scripts/peak-by-gene-pval.R $reg $m ${sample}-pval.bed

done
