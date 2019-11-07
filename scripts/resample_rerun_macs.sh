bams=$(ls /data/diana/DKS13_YNBChIPseq/OD*/*_sorted.bam)
for b in $bams; do
	
	base=$(basename $b | sed 's/_sorted.bam//' | sed 's/_YNB//' | sed 's/DKS13_//')
	echo $base
	frac=50
	newbam=/data/daniela/Projects/fixpeaks/resampled-bams/resampled_${base}.bam
	samtools view -b -s 42.$frac $b -o $newbam
	echo 'resampled bam extraction complete'
	od=$(echo $base | sed -r 's/q?[CNL][1-2]_//')
	ctrl=/data/diana/DKS13_YNBChIPseq/$od/DKS13_WCE_${od}_YNB_sorted.bam
	macs2 callpeak -t $newbam -c $ctrl -f BAM -g 18916112 -n $base -B -q 0.01 --slocal 3000 --nomodel --extsize 147 --call-summits --outdir /data/daniela/Projects/fixpeaks/MACS-resampled-bams/$base
done
