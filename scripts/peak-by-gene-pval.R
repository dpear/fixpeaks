library(GenomicRanges)
library(scales)

args         <- commandArgs(trailingOnly = TRUE)

regions      <- read.table(args[1],header=FALSE)
colnames(regions) <- c('chr','start','end','name','score','strand')
regions_granges   <- GRanges(regions)

MACS         <- read.table(args[2],header=TRUE)
MACS_granges <- GRanges(MACS)

ol <- findOverlaps(regions_granges,MACS_granges)
regions_granges@elementMetadata$log10pvals <- '-'
regions_granges@elementMetadata$log10pvals[as.numeric(ol@from)] <- MACS_granges@elementMetadata$X.log10.pvalue.[ol@to]

df <- data.frame(regions_granges)[,c(1,2,3,6,8,5)]


write.table(df,file=args[3],quote=FALSE,row.names=FALSE,col.names=FALSE)
