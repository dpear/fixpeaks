# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Rsubread")

library(edgeR)
library(Rsubread)

args <- commandArgs(trailingOnly = TRUE)
args <- c('/data/diana/DKS13_YNBChIPseq/OD1/DKS13_qN1_OD1_YNB_sorted.bam',
          '/data/diana/DKS13_YNBChIPseq/OD1/DKS13_qN2_OD1_YNB_sorted.bam',
          '/data/diana/DKS13_YNBChIPseq/OD1/DKS13_K_OD1_YNB_sorted.bam',
          'qN-edgeR-peaks-all.tsv')

saf.file  <- '~/Projects/fixpeaks/promoters.saf'
saf       <- read.table(saf.file,header=TRUE)
saf$Start[saf$Start<0] <- 0

# Link 1 &  Link 2: 2.3 Producing a table of read counts
counts <- featureCounts(file                = args[1:3],
                        annot.ext           = saf,
                        isGTFAnnotationFile = FALSE)

y <- DGEList( counts = counts$counts, group = c(1,1,2) )

# Filtering, Normalization
keep <- filterByExpr(y)
y    <- y[keep, , keep.lib.sizes=FALSE]
y    <- calcNormFactors(y)
y    <- estimateDisp(y)
e    <- exactTest(y, pair=1:2, dispersion = y$common.dispersion)
t    <- topTags(e,n= sum(keep))
t    <- t$table

t <- t[order(t$logFC),]


# View(t[t$logFC< -2 & t$PValue<.00005, ])
     
write.table(t,file=args[4],quote=FALSE,row.names = TRUE,
            col.names = c('gene_logFC','logCPM','PValue','FDR'))
print(args[4])
