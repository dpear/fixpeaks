#!/bin/bash

files=($(ls ~/Projects/fixpeaks/peak-by-gene-outputs/*))
echo 'CNAG chr start end strand' > merged
cat ${files[0]} | sort -k 4,4 | awk '{print $4, $1, $2, $3, $6}' >> merged

for f in ${files[@]}; do
	sample=$(basename $f | grep -oP 'q?[CNL][0-9]_OD[0-9]')
	header=$(echo ${sample}_rank chr start end CNAG ${sample}_height strand)
	echo $header > tmp2
	cat $f | sort -k 5,5 -r -n | awk '{printf("%d %s\n", NR, $0)}' | head -n -1 >> tmp2
	cat tmp2 | awk '{print $5,$1,$6}' | sort -k 1,1 > tmp3
	join merged tmp3 > out
	cat out > merged
done
