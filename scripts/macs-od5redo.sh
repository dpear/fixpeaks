#!/bin/bash

# for some reason (q)C2 OD5 didn't run with any of the cutoffs so here we go;
# will fix in original script 
bg=/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph
macs2 bdgpeakcall -c 3 -i $bg -o bdg3_bed/C2_OD5_3bdg.bed
macs2 bdgpeakcall -c 4 -i $bg -o bdg3_bed/C2_OD5_4bdg.bed
macs2 bdgpeakcall -c 5 -i $bg -o bdg3_bed/C2_OD5_bdg.bed
bg=/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph
macs2 bdgpeakcall -c 3 -i $bg -o bdg3_bed/qC2_OD5_3bdg.bed
macs2 bdgpeakcall -c 4 -i $bg -o bdg3_bed/qC2_OD5_4bdg.bed
macs2 bdgpeakcall -c 5 -i $bg -o bdg3_bed/qC2_OD5_bdg.bed

wc -l bdg*_bed/C2_OD5*
wc -l bdg*_bed/qC2_OD5*
