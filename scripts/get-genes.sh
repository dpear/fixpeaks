#/bin/bash

cat pval-peak-by-gene/L1_OD1-pval.bed | sort -k 5,5 -nr | awk '{if ($5!="-"){ print $4,"Liv3",$5}}' > Liv3
