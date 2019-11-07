#!/bin/bash

./scripts/cluster-peaks-by-genes.sh
# now all files should be in a folder called enriched

for f in $(ls enriched/*); do
	sample=$(basename $f | grep -Po 'q?[CLN][0-9]_OD[0-9]|WCE_OD[0-9]|K_OD[0-9]')
	awk '{print $4,FILENAME}' $f | tail -n +3 | sed 's/enriched\///' | sed -r 's/_enriched\.bed//' > raw/${sample}_raw
done
