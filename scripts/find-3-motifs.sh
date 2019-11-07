#!/bin/bash
# input: 3 motifs, cushion size
# output: bed file and fasta file with sequences containing all three motifs within a certain <cushion> of base pairs

genome=/data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta

scanMotifGenomeWide.pl C-motif.motif $genome -bed -keepAll > C-motif-genomewide.bed
cat C-motif-genomewide.bed | awk '{print $1 "\t" $10-50 "\t" $11+50 "\t" $1":"$10-50"-"$11+50 "\t" "." "\t" $14}' > C-motif-cushioned.bed

bedtools getfasta -fo C-motif-cushioned.fasta -name -fi /data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta -bed C-motif-cushioned.bed 


cat C-motif-cushioned.fasta | grep AGGAAC | grep -P 'TT[AGTC]{4}TT' | grep CAGGG

f1='AGGAAC'
r1='GTTCCT'
f2='CAGGG'
r2='CCCTG'
f3='TT[AGTC]{4}TT'
r3='AA[AGTC]{4}AA'

cat C-motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $f2 | grep -P -B 1 $f3 > all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $f2 | grep -P -B 1 $r3 >> all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $r2 | grep -P -B 1 $f3 >> all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $r2 | grep -P -B 1 $r3 >> all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $f2 | grep -P -B 1 $f3 >> all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $f2 | grep -P -B 1 $r3 >> all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $r2 | grep -P -B 1 $f3 >> all-motifs.fasta
cat C-motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $r2 | grep -P -B 1 $r3 >> all-motifs.fasta

cat all-motifs.fasta | grep -vP '\-\-' > temp
mv temp all-motifs.fasta

awk '/^>/{f=!d[$1];d[$1]=1}f' all-motifs.fasta > deduped-all-motifs.fasta
cat deduped-all-motifs.fasta | grep -P '\>' | sed 's|[>]||g' > sequence-headers.txt
cat C-motif-cushioned.bed | grep -f sequence-headers.txt > all-motifs-50.bed


####################### N MOTIF #######################

genome=/data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta

scanMotifGenomeWide.pl N-motif.motif $genome -bed -keepAll > motif-genomewide.bed
cat motif-genomewide.bed | awk '{print $1 "\t" $10-25 "\t" $11+25 "\t" $1":"$10-25"-"$11+25 "\t" "." "\t" $14}' > motif-cushioned.bed

bedtools getfasta -bed motif-cushioned.bed -fo motif-cushioned.fasta -name -fi /data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta

f1='.'
r1='.'
f2='AGGG'
r2='CCCT'
f3='TT[AGTC]{4}TT'
r3='AA[AGTC]{4}AA'

cat motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $f2 | grep -P -B 1 $f3 > all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $f2 | grep -P -B 1 $r3 >> all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $r2 | grep -P -B 1 $f3 >> all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $f1 | grep -P -B 1 $r2 | grep -P -B 1 $r3 >> all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $f2 | grep -P -B 1 $f3 >> all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $f2 | grep -P -B 1 $r3 >> all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $r2 | grep -P -B 1 $f3 >> all-motifs.fasta
cat motif-cushioned.fasta | grep -P -B 1 $r1 | grep -P -B 1 $r2 | grep -P -B 1 $r3 >> all-motifs.fasta

cat all-motifs.fasta | grep -vP '\-\-' > temp
mv temp all-motifs.fasta

awk '/^>/{f=!d[$1];d[$1]=1}f' all-motifs.fasta > deduped-all-motifs.fasta
cat deduped-all-motifs.fasta | grep -P '\>' | sed 's|[>]||g' > sequence-headers.txt
cat motif-cushioned.bed | grep -f sequence-headers.txt > all-motifs-short.bed


