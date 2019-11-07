#!/bin/bash

macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/C1_OD1_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/C1_OD5_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/C2_OD1_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/C2_OD5_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/qC1_OD1_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC1_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/qC1_OD5_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/qC2_OD1_bdg.bed
macs2 bdgpeakcall	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg/qC2_OD5_bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/C1_OD1_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/C1_OD5_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/C2_OD1_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/C2_OD5_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/qC1_OD1_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC1_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/qC1_OD5_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/qC2_OD1_3bdg.bed
macs2 bdgpeakcall -c 3	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg3/qC2_OD5_3bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/C1_OD1_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/C1_OD5_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/C2_OD1_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/C2_OD5_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/qC1_OD1_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC1_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/qC1_OD5_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/qC2_OD1_4bdg.bed
macs2 bdgpeakcall -c 4	-i	/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qC2_OD5_YNB_sorted_norm_50bp_smooth.bedgraph	-o	../bed_bdg4/qC2_OD5_4bdg.bed
