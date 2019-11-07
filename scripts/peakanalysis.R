############ Installation etc. ############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("Rsamtools")
# https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html

library(ChIPpeakAnno)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)

###########################################

time1 <- read.delim("/data/diana/DKS11and12_ChIPseqTimecourse/MACSanalysis/L_logA_sorted/L_logA_sorted_peaks.xls",skip=28,header=TRUE)
time2 <- read.delim("/data/diana/DKS11and12_ChIPseqTimecourse/MACSanalysis/L_OD5A_sorted/L_OD5A_sorted_peaks.xls",skip=28,header=TRUE)
time3 <- read.delim("/data/diana/DKS11and12_ChIPseqTimecourse/MACSanalysis/L_OD10A_sorted/L_OD10A_sorted_peaks.xls",skip=28,header=TRUE)


time1range <- rangeformat(time1)
time2range <- rangeformat(time2)
time3range <- rangeformat(time3)

View(time3)

DATA_DIRECTORY <- "/data/diana/DKS11and12_ChIPseqTimecourse/MACSanalysis"

xlss <- read.delim(DATA_DIRECTORY)

# Get data from macs output, which i am assuming is somewhat correct
xlss       <- list.files("MACSoutput/",full.names=TRUE)
readxls    <- function(x) read.delim(x,header=TRUE,colClasses = c("character",rep("integer",4),rep("double",4),"character"))
xls        <- lapply(xlss,readxls)
names(xls) <- list.files("MACSoutput/")

##### functions ########
rangeformat       <- function(A){
  repA  <- toGRanges(A,format="MACS",header=TRUE)
  repA@ranges@NAMES <- paste("",seq(1:length(repA)))
  return(repA)
} # A: .xls dataframes read in already. returns granges object with metadata.
addEnrichment     <- function(fullFrame){
  qry <- fullFrame
  sbj <- reduce(fullFrame)
  
  `%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
  
  hts <- findOverlaps(qry, sbj)
  sbj$fold_enrichment <- rep(0,length(sbj))
  sbj$hits <- rep(0,length(sbj))
  
  qry$hits <- rep(1,length(qry))
  sbj$fold_enrichment[hts@to] %+=% qry$fold_enrichment[hts@from]
  sbj$hits[hts@to]     %+=% qry$hits[hts@from]
  sbj$fold_enrichment <- sbj$fold_enrichment / sbj$hits
  return(sbj) 
}  # fullFrame: output from rangeformat. returns overlappingpeaks object with metadata attached cleverly.
compareReplicates <- function(A,B){
  # convert into reduced granges with no overlaps
  repA  <- rangeformat(A)
  repB  <- rangeformat(B)
  repA <- addEnrichment(repA)
  repB <- addEnrichment(repB)
  
  # look at overlapping metadata (trying to encorporate metadata)
  ol_AB <- findOverlapsOfPeaks(repA,repB,maxgap = 1000)
  makeVennDiagram(ol_AB,fill=c("steelblue", "blue3"))
  return(ol_AB)
} # A,B: .xls dataframes read in already. returns overlaps of replicates a and b without metadata (boo)

########################

# looking at L- last time period
L  <- compareReplicates(xls$L03A.xls,xls$L03B.xls)
L_A <- rangeformat(xls$L03A.xls)
L_B <- rangeformat(xls$L03B.xls)
qL <- compareReplicates(xls$qL03A.xls,xls$qL03B.xls)
qL_A <- rangeformat(xls$qL03A.xls)
qL_B <- rangeformat(xls$qL03B.xls)

# check to make sure the enrichment of merged peaks isnt different
hist(L$peaksInMergedPeaks$fold_enrichment)
hist(L_A$fold_enrichment)
hist(L_B$fold_enrichment)

hist(qL$peaksInMergedPeaks$fold_enrichment)
hist(qL_A$fold_enrichment)
hist(qL_B$fold_enrichment)
# good to go

L$peaksInMergedPeaks$origin <- rep('L',length(L$peaksInMergedPeaks))
qL$peaksInMergedPeaks$origin <- rep('qL',length(qL$peaksInMergedPeaks))
overlap_L <- findOverlapsOfPeaks(L$peaksInMergedPeaks,qL$peaksInMergedPeaks)


p <- ggplot(data=data.frame(overlap_L$uniquePeaks))+
  geom_point(aes(x=start,col=fold_enrichment,y=factor(origin),alpha=.4))+
  xlim(0,500000)
p

# Error: incorrect number of dimensions
L_unique <- overlap_L$uniquePeaks$fold_enrichment[overlap_L$uniquePeaks$origin=='L',]
qL_unique <- overlap_L$uniquePeaks$fold_enrichment[overlap_L$uniquePeaks$origin=='qL',]

p <- ggplot() + 
  geom_density(data=data.frame(L$peaksInMergedPeaks),aes(x=fold_enrichment,alpha=.4,col='blue')) +
  geom_density(data=data.frame(qL$peaksInMergedPeaks),aes(x=fold_enrichment,alpha=.4,col='black'))
p



qL <- compareReplicates(xls$qL03B.xls,xls$qL03B.xls)
N  <- compareReplicates(xls$N03A.xls,xls$N03B.xls)
qN <- compareReplicates(xls$qN03A.xls,xls$qN03B.xls)

ol_L_qL <- findOverlapsOfPeaks(L$mergedPeaks,qL$mergedPeaks,maxgap = 1000)
ol_N_qN <- findOverlapsOfPeaks(N$mergedPeaks,qN$mergedPeaks,maxgap = 1000)
makeVennDiagram(ol_N_qN)

### look at overall enrichment of peaks of all samples
par(mfrow=(c(1,1)))
hist(1)
getEnrichment      <- function(x) return(x$fold_enrichment)
enrichments        <- lapply(xls,getEnrichment)
names(enrichments) <- list.files("MACSoutput/")
par(mfrow=c(4,6))
hists <- lapply(enrichments,hist)
for (i in 1:length(hists)) hists[[i]]$xname=names(xls)[i]
hist(1)
lapply(hists,plot)
hist(1)


