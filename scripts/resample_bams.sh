#!/bin/bash

# resamples from files listed in bams, extracts <min_reads> from each bamfile and outputs

bams=$(ls /data/diana/DKS13_YNBChIPseq/OD*/*_sorted.bam)

min_reads=9942317

for b in $bams; do

        base=$(basename $b | sed 's/_sorted.bam//' | sed 's/_YNB//' | sed 's/DKS13_//')
        echo $base
	reads=$(samtools view -c $b)
        frac=$(echo $(( 100 * $min_reads / $reads ))) # cheating in unix floating point division
        newbam=/data/daniela/Projects/fixpeaks/resampled-bams/resampled_${base}.bam
        samtools view -b -s 42.$frac $b -o $newbam
        echo 'resampled bam extraction complete'
       
done
