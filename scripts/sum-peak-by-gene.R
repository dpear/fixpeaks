library(GenomicRanges)
library(dplyr)
# args <- c('/data/diana/DKS13_YNBChIPseq/Bedgraphs/DKS13_C1_OD1_YNB_sorted_norm_50bp_smooth.bedgraph','~/Projects/fixpeaks/promoters.txt','~/Projects/fixpeaks/C1_OD1_sums.bed')

# read files and load to Granges
args            <- commandArgs(trailingOnly = TRUE)
bedgraph        <- read.table(args[1],header=FALSE,sep="\t")
genes           <- read.table(args[2],header=FALSE,sep="\t")
names(bedgraph) <- c('chr','start','end','height')
names(genes)    <- c('chr','start','end','gene','score','strand')
bedgraph        <- GRanges(bedgraph)
genes           <- GRanges(genes)

# assign genes to bedgraph file
ol <- findOverlaps(genes, bedgraph)
bedgraph@elementMetadata$gene <- '-'
bedgraph@elementMetadata$gene[ol@to] <- as.character(genes@elementMetadata$gene[ol@from]) 

# get sums of each gene
bedgraph    <- data.frame(bedgraph)
sums        <- aggregate(bedgraph$height, by=list(bedgraph$gene), sum)
names(sums) <- c('gene', 'score')

# join the sum and genes data
both <- inner_join(sums, data.frame(genes), by=c('gene'))
both <- both[ ,c(3,4,5,1,2,7)]

write.table(both,args[3],row.names = FALSE,col.names = FALSE,quote=FALSE)
