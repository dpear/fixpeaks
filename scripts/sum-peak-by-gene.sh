#!/bin/bash

bgs=(/data/diana/DKS14_TFmutChIPseq/DKS14_round2_SCHM28/DKS14_cL1-b_YNB_sorted_norm_50bp_smooth.bedgraph /data/diana/DKS14_TFmutChIPseq/DKS14_round2_SCHM28/DKS14_cL2-b_YNB_sorted_norm_50bp_smooth.bedgraph)
genes=/home/daniela/Projects/fixpeaks/promoters.txt
for bg in ${bgs[@]}; do
	sample=$(echo $bg | grep -Po 'q?c?n?l?[CLN][0-9]-b|WCE[0-9]')
	Rscript scripts/sum-peak-by-gene.R $bg $genes ${sample}_sums.bed
	echo $sample
done
	
