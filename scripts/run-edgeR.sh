#!/bin/bash

DATA=/data/diana/DKS13_YNBChIPseq/OD1/
untagged=/data/diana/DKS13_YNBChIPseq/OD1//DKS13_K_OD1_YNB_sorted.bam
Rscript run-edgeR.R $DATA/DKS13_C1_OD1_YNB_sorted.bam $DATA/DKS13_C2_OD1_YNB_sorted.bam $untagged C-edgeR-peaks.tsv
Rscript run-edgeR.R $DATA/DKS13_N1_OD1_YNB_sorted.bam $DATA/DKS13_N2_OD1_YNB_sorted.bam $untagged N-edgeR-peaks.tsv
Rscript run-edgeR.R $DATA/DKS13_L1_OD1_YNB_sorted.bam $DATA/DKS13_L2_OD1_YNB_sorted.bam $untagged L-edgeR-peaks.tsv
Rscript run-edgeR.R $DATA/DKS13_qC1_OD1_YNB_sorted.bam $DATA/DKS13_qC2_OD1_YNB_sorted.bam $untagged qC-edgeR-peaks.tsv
Rscript run-edgeR.R $DATA/DKS13_qN1_OD1_YNB_sorted.bam $DATA/DKS13_qN2_OD1_YNB_sorted.bam $untagged qN-edgeR-peaks.tsv
Rscript run-edgeR.R $DATA/DKS13_qL1_OD1_YNB_sorted.bam $DATA/DKS13_qL2_OD1_YNB_sorted.bam $untagged qL-edgeR-peaks.tsv
