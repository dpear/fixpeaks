# Link 1: # https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts
# Link 2: # https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")

library(edgeR)
library(Rsubread)
args <- commandArgs(trailingOnly = TRUE)


saf.file  <- '~/Projects/fixpeaks/promoters.saf'

bam.files <- args[1:4]

saf       <- read.table(saf.file,header=TRUE)
saf$Start[saf$Start<0] <- 0

# Link 1 &  Link 2: 2.3 Producing a table of read counts
counts <- featureCounts(file                = bam.files,
                        annot.ext           = saf,
                        isGTFAnnotationFile = FALSE)

# Link 2: 
# skipping normalization to untagged:
# "any technical factor that is unrelatedto the experimental conditions 
#  should cancel out of any differential expression analysis"

y <- DGEList(counts = counts$counts,
             group = c(1,1,2,2))


#y <- DGEList(counts = counts$counts,
             group  = c(1,1, 2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,
                        9,9,10,11,11,12,12,13,13,14,15,15,16))
#y <- DGEList(counts = counts$counts,
             group  = c(1,1, 2, 3,3, 4,4, 5,5, 6,6, 7,7, 8))
#y$samples

# 2.6 Filtering
# filter out genes with too low exptression
keep <- filterByExpr(y)
y    <- y[keep, , keep.lib.sizes=FALSE]
#sum(keep) # number of genes kept

# Link 2: 2.7.3
y <- calcNormFactors(y)
#y$samples
# normalization factor below one = small number of high count genes are dominating
# causing counts for other genes to be lower than usual given library size

y <- estimateDisp(y)

e <- exactTest(y, pair=1:2, dispersion = y$common.dispersion)
t <- topTags(e,n= sum(keep))


#e <- exactTest(y,pair=c(1,5))
#t <- topTags(e,n= 6871)
t <- t$table
write.table(t,file=args[5],quote=FALSE,row.names = TRUE,
            col.names = c('gene_logFC','logCPM','PValue','FDR'))

#e <- exactTest(y,pair=c(3,6))
#t <- topTags(e,n= 6871)
#t <- t$table
#write.table(t,file='~/Projects/fixpeaks/L_OD1-edgeR.tsv',quote=FALSE,row.names = TRUE,
#            col.names = c('gene_logFC','logCPM','PValue','FDR'))

#e <- exactTest(y,pair=c(4,7))
#t <- topTags(e,n= 6871)
#t <- t$table
#write.table(t,file='~/Projects/fixpeaks/N_OD1-edgeR.tsv',quote=FALSE,row.names = TRUE,
#            col.names = c('gene_logFC','logCPM','PValue','FDR'))





# 
# # ~~~~ graveyard ~~~~~~
# y <- estimateDisp(y)
# # y <- estimateCommonDisp(y)
# # y <- estimateTagwiseDisp(y)
# 
# groupings   <- list(c('DKS13.C1.OD1.YNB.sorted.bam','DKS13.qC1.OD5.YNB.sorted.bam'),
#                     c('DKS13.C2.OD1.YNB.sorted.bam','DKS13.qC2.OD5.YNB.sorted.bam'),
#                     c('DKS13.L1.OD1.YNB.sorted.bam','DKS13.qL1.OD5.YNB.sorted.bam'),
#                     c('DKS13.L2.OD1.YNB.sorted.bam','DKS13.qL2.OD5.YNB.sorted.bam'),
#                     c('DKS13.N1.OD1.YNB.sorted.bam','DKS13.qN1.OD5.YNB.sorted.bam'),
#                     c('DKS13.N2.OD1.YNB.sorted.bam','DKS13.qN2.OD5.YNB.sorted.bam'),
#                     
#                     c('DKS13.C1.OD5.YNB.sorted.bam','DKS13.qC1.OD5.YNB.sorted.bam'),
#                     c('DKS13.C2.OD5.YNB.sorted.bam','DKS13.qC2.OD5.YNB.sorted.bam'),
#                     c('DKS13.L1.OD5.YNB.sorted.bam','DKS13.qL1.OD5.YNB.sorted.bam'),
#                   
#                     c('DKS13.N1.OD5.YNB.sorted.bam','DKS13.qN1.OD5.YNB.sorted.bam'),
#                     c('DKS13.N2.OD5.YNB.sorted.bam','DKS13.qN2.OD5.YNB.sorted.bam'))
#   
#   
#   
#   
# comparisons <- list()
# 
# n <- length(groupings)
# for (i in 1:n){
#  comparisons[[i]]       <- exactTest(y,pair=groupings[[i]])
#  comparisons[[i]]$table <- comparisons[[i]]$table[order(comparisons[[i]]$table$PValue),]
# }
# 
# # Link 2: 2.11 What to do if you have no replicates
# # Just for Liv3 OD5
# comparisons[[n+1]]       <- exactTest(y,pair=c(11,14),dispersion = y$common.dispersion)
# comparisons[[n+1]]$table <- comparisons[[n+1]]$table[order(comparisons[[n+1]]$table$PValue),]
# 
# n <- nrow(comparisons[[1]]$table)
# t <- list()
# for (i in 1:length(comparisons)){
#   t      <- topTags(comparisons[[i]],n=n)
#   t[[i]] <- t[t$table$FDR<.0001,]
# }
# 
# et <- exactTest(y,pair=c(1,5))
# topTags(et)
# et
# 
# 
# # Look at 3 tf in WT
# # does not work only returns peakless genes
# # one-way ANOVA test
# design <- model.matrix(~group, data=y$samples)
# design
# 
# fit <- glmQLFit(y, design)
# qlf <- glmQLFTest(fit, coef=3:4)
# topTags(qlf)