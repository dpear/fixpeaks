#!/bin/bash

python3 locate-motif.py /data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta ../CCCT-hits.bed 10
bedtools getfasta -fo ../CCCT-hits.fasta -fi /data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta -bed ../CCCT-hits.bed
cat ../CCCT-hits.fasta | grep -P -B 1 'TT[AGTC]{4}TT|AA[ACTC]{4}AA' | grep -P '\>' | sed 's/>//' | sed -r 's/\:/\t/' | sed -r 's/-/\t/' > ../Nrg1_Liv3_hits.bed
