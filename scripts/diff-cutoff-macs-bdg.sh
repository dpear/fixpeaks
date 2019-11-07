#!/bin/bash

files=$(ls /data/diana/DKS13_YNBChIPseq/Bedgraphs/*norm_50bp_smooth.bedgraph )

for f in $files; do
	echo $(basename $f)
	avg_ht=$(cat $f | awk '{ total += $4 } END { print total/NR }')
	echo $avg_ht





done;
