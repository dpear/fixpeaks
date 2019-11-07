#!/bin/bash

# usage: ./sicer-analysis.sh 

# convert bam to bed
# comment out once this has been done once

bams=($(ls /data/diana/DKS13_YNBChIPseq/*/*_sorted.bam))
windowsize=200

for bam in ${bams[@]}; do
	sample=$(echo $bam | grep -oP 'q?[CLN][0-9]_OD[0-9]|K_OD[0-9]|WCE_OD[0-9]')
	bedtools bamtobed -i $bam > ${sample}_frombam.bed
done

# sh /home/jordan/PROGRAMS/sicer/SICER-df.sh 
