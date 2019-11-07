#!/bin/bash
f2=$(ls /data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_WCE*sorted_norm_50bp_smooth.bedgraph)
f1=$(ls /data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qN2*sorted_norm_50bp_smooth.bedgraph)
bgs=$(echo $f1 $f2)
echo $bgs

for bg in $bgs; do

	sample=$(echo $bg | grep -Po 'q?[CLN][0-9]_OD[0-9]|WCE_OD[0-9]')
	echo $sample

	# grab window upstream (depends on strand)
	cat $genes | grep '-' | awk -v size=$window '{print $4 "\t" $3 "\t" $3+size "\t" $1 "\t" "score" "\t" $5}' > tmp
	cat $genes | grep '+' | awk -v size=$window '{print $4 "\t" $2-size "\t" $2 "\t" $1 "\t" "score" "\t" $5}' >> tmp
	echo 'gene ranges edited'

	# get the max and min position of each chr recorded by the bedgraph
	chrs=$(cat $genes | awk '{print $4}' | uniq | grep chr)
	echo '' > tmp2
	for chr in $chrs; do 
		min=$(cat $bg | grep -P "${chr}\t" | head -1 | awk '{print $2}') 
		max=$(cat $bg | grep -P "${chr}\t" | tail -1 | awk '{print $3}') 
		cat tmp | grep -P "${chr}\t" | awk -v min=$min -v max=$max '{print $0 "\t" min "\t" max}'
	done >> tmp2
	echo 'chr max and mins extracted'

	# correct any out-of-bounds regions created by
	# looking at 1kb upstream of genes
	cat tmp2 | awk '\
	function min(a,b) {return a < b ? a: b}
	function max(a,b) {return a > b ? a: b}
	{print $1 "\t" max($2,$7) "\t" min($3,$8) "\t" $4 "\t" $5 "\t" $6}' > tmp3
	echo 'promoter region file created'

	# delete any rows with starts>ends
	cat tmp3 | awk '{ if($2<$3) {print $0}}' > tmp_f; cat tmp_f > tmp3

	echo '' > ${sample}_final.bed
	for chr in $chrs; do
		cat $bg | grep -P "${chr}\t" > bg-tmp # bg heights for this chr
		cat tmp3 | grep -P "${chr}\t" > p-tmp # p-tmp = genes on this chr
		echo $chr 1 # reference while running
		cat p-tmp | while read line; do 
			min=$(echo $line | awk '{print $2}')
			max=$(echo $line | awk '{print $3}')
			score=$(cat bg-tmp | awk -v min=$min -v max=$max '{if($2>=min && $3<=max) {print $0}}' | sort -n -k 4,4 -r | head -1 | awk '{print $4}')
			echo $line | awk -v score=$score '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" score "\t" $6}' >> ${sample}_final.bed
			cat final | tail -1
		done
		echo $chr 2 # reference while running
	done
done



rm tmp
rm tmp2
rm tmp3
rm bg-tmp
rm p-tmp
