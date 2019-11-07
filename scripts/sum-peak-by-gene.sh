#!/bin/bash

bgs=$(ls /data/diana/DKS13_YNBChIPseq/Bedgraphs/*sorted_norm_50bp_smooth.bedgraph| grep -v normWCE| grep -v WCE)
bgs=$(ls
genes=/home/daniela/Projects/fixpeaks/promoters.txt

for bg in ${bgs[@]}; do
	sample=$(echo $bg | grep -Po 'q?[CLN][0-9]_OD[0-9]|WCE_OD[0-9]|K_OD[0-9]')
	Rscript scripts/sum-peak-by-gene.R $bg $genes ${sample}_sums.bed
done
	
