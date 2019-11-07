#!/bin/bash

max=/data/daniela/Projects/fixpeaks/resampled-bams/resampled_qL1_OD1.bam
max_bg=/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qL1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph
min=/data/daniela/Projects/fixpeaks/resampled-bams/resampled_qN2_OD1.bam
min_bg=/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_qN2_OD1_YNB_sorted_norm_50bp_smooth.bedgraph
ctrl=/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_WCE_OD1_YNB_sorted_norm_50bp_smooth.bedgraph


macs2	callpeak	-t	$max	-c	$ctrl	--slocal	1000			-B	-f	BAM	-g	18899812	--call-summits	--outdir	1	-n	1	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	1000			-B	-f	BAM	-g	18899812		--outdir	2	-n	2	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	1000	--max-gap	0	-B	-f	BAM	-g	18899812	--call-summits	--outdir	3	-n	3	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	1000	--max-gap	0	-B	-f	BAM	-g	18899812		--outdir	4	-n	4	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	2000			-B	-f	BAM	-g	18899812	--call-summits	--outdir	5	-n	5	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	2000			-B	-f	BAM	-g	18899812		--outdir	6	-n	6	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	2000	--max-gap	0	-B	-f	BAM	-g	18899812	--call-summits	--outdir	7	-n	7	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	2000	--max-gap	0	-B	-f	BAM	-g	18899812		--outdir	8	-n	8	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	3000			-B	-f	BAM	-g	18899812	--call-summits	--outdir	9	-n	9	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	3000			-B	-f	BAM	-g	18899812		--outdir	10	-n	10	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	3000	--max-gap	0	-B	-f	BAM	-g	18899812	--call-summits	--outdir	11	-n	11	--nomodel	--extsize	75
macs2	callpeak	-t	$max	-c	$ctrl	--slocal	3000	--max-gap	0	-B	-f	BAM	-g	18899812		--outdir	12	-n	12	--nomodel	--extsize	75
macs2	bdgpeakcall	-i	$max_bg													--outdir	13		13			
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	1000			-B	-f	BAM	-g	18899812	--call-summits	--outdir	14	-n	14	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	1000			-B	-f	BAM	-g	18899812		--outdir	15	-n	15	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	1000	--max-gap	0	-B	-f	BAM	-g	18899812	--call-summits	--outdir	16	-n	16	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	1000	--max-gap	0	-B	-f	BAM	-g	18899812		--outdir	17	-n	17	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	2000			-B	-f	BAM	-g	18899812	--call-summits	--outdir	18	-n	18	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	2000			-B	-f	BAM	-g	18899812		--outdir	19	-n	19	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	2000	--max-gap	0	-B	-f	BAM	-g	18899812	--call-summits	--outdir	20	-n	20	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	2000	--max-gap	0	-B	-f	BAM	-g	18899812		--outdir	21	-n	21	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	3000			-B	-f	BAM	-g	18899812	--call-summits	--outdir	22	-n	22	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	3000			-B	-f	BAM	-g	18899812		--outdir	23	-n	23	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	3000	--max-gap	0	-B	-f	BAM	-g	18899812	--call-summits	--outdir	24	-n	24	--nomodel	--extsize	75
macs2	callpeak	-t	$min	-c	$ctrl	--slocal	3000	--max-gap	0	-B	-f	BAM	-g	18899812		--outdir	25	-n	25	--nomodel	--extsize	75
macs2	bdgpeakcall	-i	$min_bg													--outdir	26					
