#!/bin/bash

bams=($(ls /data/diana/DKS13_YNBChIPseq/OD*/*_sorted.bam | grep -v K_ | grep -v WCE))

for bam in ${bams[@]}; do
	sample=$(basename $bam | grep -Po 'q?[CLN][0-9]_OD[0-9]|WCE_OD[0-9]|K_OD[0-9]')
	OD=$(basename $bam | grep -Po '_OD[0-9]'| sed 's/_//')
	ctrl=/data/diana/DKS13_YNBChIPseq/${OD}/DKS13_K_${OD}_YNB_sorted.bam
	extsize=$(cat fragment-sizes.txt | grep -P "^$sample" | awk '{print $2}')
	echo $sample
	echo $OD
	echo $ctrl
	macs2 callpeak -t $bam -g 18899812 -c $ctrl --nomodel --extsize $extsize --outdir ${sample}_K_control
done
