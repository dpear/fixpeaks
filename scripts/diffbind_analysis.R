# DiffBind documentation can be found at http://169.230.6.100:8787/help/library/DiffBind/doc/DiffBind.pdf
# author: Daniela Perry danielaperry2015@gmail.com
# installation
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install('grimbough/Rhtslib')
BiocManager::install("Rsamtools")
BiocManager::install("DiffBind")
BiocManager::install("biomaRt")

library(Rsamtools)
library(DiffBind)
library(GenomicRanges)
library(ChIPpeakAnno)
library(scales)
library(dplyr)
library(DESeq2)

samples <- read.csv("~/Projects/fixpeaks/samplesheets/bdg-sheet-all.csv", sep="\t")
peaks.d <- dba(sampleSheet = samples)
peaks   <- dba.count(peaks.d)
peaks   <- dba.contrast(peaks, categories = DBA_CONDITION, minMembers = 2)
peaks   <- dba.analyze(peaks)
peaks.DB<- dba.report(peaks) 
peaks.DB<- dba.report(peaks,th=.01) 

peaks.rep <- dba.report(peaks,bCalled=TRUE,th=1)
peaks.rep <- peaks.rep[order(-peaks.rep$Fold^2),]
ol        <- findOverlaps(peaks.rep,genes)
peaks.rep@elementMetadata$genes <- '-'
peaks.rep@elementMetadata$genes[as.numeric(ol@from)] <- genes@elementMetadata$gene[ol@to]
peaks.rep[ol@from]@elementMetadata$genes <-  as.character(genes@elementMetadata$gene[ol@to])
peaks.rep@elementMetadata$genes <- as.factor(peaks.rep@elementMetadata$genes)
peaks.rep

df      = data.frame(flank(genes,1000))
outname = "~/Projects/fixpeaks/upstream-regions.tsv"
write.table(df,file=outname,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")



# plotting
plot(peaks.d)
dba.plotPCA(peaks)
dba.plotMA(peaks)
dba.plotVolcano(peaks)
dba.plotBox(peaks)
dba.plotHeatmap(peaks, score=DBA_SCORE_RPKM_FOLD,
                colScheme = 'RdPu') # different than above, see docs
dba.plotHeatmap(peaks, contrast=1, correlations=FALSE,th=.1,
                colScheme='RdPu',main="") # binding affinity for differentially bound sites
dba.plotVenn(peaks,1:4,label1="WT (rep A)",label2="WT (rep B)",
             label3="Mutant (repA)", label4="Mutant (rep B)")


# occupancy analysis
peaks.rep <- dba.report(peaks,bCalled=TRUE,th=1)
only.wt   <- peaks.rep$Called1>0 & peaks.rep$Called2==0
sum(only.wt)
only.mut  <- peaks.rep$Called1==0 & peaks.rep$Called2>0
sum(only.mut)
both      <- peaks.rep$Called1>0 & peaks.rep$Called2>0
sum(both)


#### COLUMNS IN DB~ from http://169.230.6.100:8787/help/library/DiffBind/doc/DiffBind.pdf
# The value columns show the mean read concentration over all the samples 
# (the default calculation uses log2 normalized ChIP read counts 
# with control read counts subtracted) and the mean concentration 
# over the first (Resistant) group and second (Responsive) group. 
# The Fold column shows the difference in mean concentrations 
# between the two groups (Conc_Resistant - Conc_Responsive), 
# with a positive value indicating increased binding affinity in the Resistant group 
# and a negative value indicating increased binding affinity in the Responsive group.
# The final two columns give confidence measures for identifying these sites 
# as differentially bound, with a raw p-value and a multiple testing corrected FDR 
# in the final column

# view genes from modified gff file and create granges object
genes <- read.csv("/data/daniela/Projects/fixpeaks/files/h99-cneoformans-noncoding-removed.gff",sep="\t")
genes <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
genes

# change the start (+ strand) or end (- strand) to 500 bp upstream of the TSS
end(genes[genes@strand=='-']) <- end(genes[genes@strand=='-']) + 500
start(genes[genes@strand=='+']) <- start(genes[genes@strand=='+']) - 500
genes

# assign gene names to peaks.DB (see above) 
ol <- findOverlaps(peaks.DB,genes)
peaks.DB@elementMetadata$genes <- '-'
peaks.DB@elementMetadata$genes[as.numeric(ol@from)] <- genes@elementMetadata$gene[ol@to]
peaks.DB[ol@from]@elementMetadata$genes <-  as.character(genes@elementMetadata$gene[ol@to])
#peaks.rep@elementMetadata$genes <- as.factor(peaks.rep@elementMetadata$genes)
peaks.DB

# write the significant peaks to a bed file
s <- data.frame(peaks.DB)
s <- data.frame(s$seqnames,s$start,s$end,s$genes,s$Fold)
names(s) <- c('chrom','chromStart','chromEnd','name','score')
s$score  <- floor(10*((s$score)^2)) # score is just for visualization in igv
write.csv(s,file="Nfinal_analysis.bed",quote = FALSE,row.names = FALSE)
range(s$score)

# get genes list with lowest p-value
peaks.rep <- dba.report(peaks,bCalled=TRUE,th=1)
peaks.rep <- peaks.rep[order(-peaks.rep$Fold^2),]
ol        <- findOverlaps(peaks.rep,genes)

peaks.rep@elementMetadata$genes <- '-'
peaks.rep@elementMetadata$genes[as.numeric(ol@from)] <- genes@elementMetadata$gene[ol@to]
peaks.rep[ol@from]@elementMetadata$genes <-  as.character(genes@elementMetadata$gene[ol@to])
peaks.rep@elementMetadata$genes <- as.factor(peaks.rep@elementMetadata$genes)

peaks.a <- aggregate(FDR ~ genes, peaks.rep, min)
peaks.m <- merge(peaks.a,peaks.rep, by.x = c("genes","FDR"), by.y=c("genes","FDR"))
peaks.m <- peaks.m[order(-peaks.m$Fold^2),]
#peaks.m$wt <- peaks.m$Called1>0
#peaks.m$mut <- peaks.m$Called2>0
#peaks.m$unique <- (peaks.m$wt==TRUE & peaks.m$mut==FALSE) | (peaks.m$wt==FALSE & peaks.m$mut==TRUE)
#peaks.m[peaks.m$unique==TRUE | peaks.m$FDR<.05,]
head(peaks.m)

names(peaks.m) <- c(   "genes_c1",
                         "FDR_c1",    
                    "seqnames_c1", 
                       "start_c1",   
                         "end_c1",    
                       "width_c1", 
                      "strand_c1",  
                        "Conc_c1",   
                     "Conc_wt_c1", 
                    "Conc_mut_c1",
                        "Fold_c1", 
                     "p.value_c1",  
                     "Called1_c1", 
                     "Called2_c1"  ) 
peaks.m$in.c1=TRUE
head(peaks.m)

write.csv(peaks.m,file="C-OD1-sheet.csv",quote = FALSE,row.names = FALSE)

c <- read.csv("/data/daniela/Projects/fixpeaks/files/C-OD5_sheet.csv")
n <- read.csv("/data/daniela/Projects/fixpeaks/files/N-OD5_sheet.csv")
#l <- read.csv("/data/daniela/Projects/fixpeaks/files/L-OD5_sheet.csv")

cn <- full_join(c,n,by=c("genes_c1"="genes_n1"))




# final analysis
n <- read.csv("~/N_sheet.csv")
head(n)
l <- read.csv("~/L_sheet.csv")
head(l)
join <- full_join(n,l,by=c("genes_n"="genes_l"))
join$interesting <- join$FDR_n<.05 | join$FDR_l<.05 | join$unique_l==TRUE | join$unique_n==TRUE

keep_rows <- join[join$interesting==TRUE,]
keep_rows <- na.omit(keep_rows)
inspection <- keep_rows[,c(1,2,18,19,34,35)]

n <- n[order(-n$Fold_n^2),]




# ordering playing around
peaks.rep[order(-peaks.rep$Conc/peaks.rep@ranges@width),]
# first entry is tight , second is weird





