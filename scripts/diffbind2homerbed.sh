#!/bin/bash

# convert columns from diffbind output to homer input for motif finding and fasta sequence extraction
# mainly for my reference

cat diffbind-output/C-N-L.tsv | awk '{print $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3 "\t" "." "\t" "+"}' > homer-input/C-N-L.bed
