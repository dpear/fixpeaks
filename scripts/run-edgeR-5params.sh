#!/bin/bash

DATA=/data/diana/DKS13_YNBChIPseq/OD1/
Rscript run-edgeR.R $DATA/DKS13_C1_OD1_YNB_sorted.bam $DATA/DKS13_C2_OD1_YNB_sorted.bam $DATA/DKS13_qC1_OD1_YNB_sorted.bam $DATA/DKS13_qC2_OD1_YNB_sorted.bam C-edgeR-diff.tsv
Rscript run-edgeR.R $DATA/DKS13_N1_OD1_YNB_sorted.bam $DATA/DKS13_N2_OD1_YNB_sorted.bam $DATA/DKS13_qN1_OD1_YNB_sorted.bam $DATA/DKS13_qN2_OD1_YNB_sorted.bam N-edgeR-diff.tsv
Rscript run-edgeR.R $DATA/DKS13_L1_OD1_YNB_sorted.bam $DATA/DKS13_L2_OD1_YNB_sorted.bam $DATA/DKS13_qL1_OD1_YNB_sorted.bam $DATA/DKS13_qL2_OD1_YNB_sorted.bam L-edgeR-diff.tsv
