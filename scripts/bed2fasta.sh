#!/bin/bash
N=$(ls bdg_bed/*[NCL]*)
for f in $N; do
	base=$(basename $f | sed -r 's/bdg_bed\///' | sed -r 's/_bdg\.bed//')
	bedtools getfasta -fi /data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta -bed $f -fo bdg_fasta/${base}_bdg.fasta
done	
