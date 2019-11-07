library(scales)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) { stop('Usage: Rscript cluster.R in.bed out.bed', call.=FALSE) }

df           <- read.table(args[1], header=FALSE)
colnames(df) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
df$cluster   <- kmeans(df$score,3)$cluster
#df$score     <- rescale(df$score,to=c(0,1000))

m <- tail(df$cluster[order(-df$score)])[1]
df$cluster[df$cluster!=m] <- 0

df <- df[df$cluster==0,]
df$cluster <- NULL

write.table(df, file=args[2], row.names=FALSE, quote=FALSE, col.names = FALSE)

# chr start end CNAG height strand (.bed columns)