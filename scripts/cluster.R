library(scales)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) { stop('Usage: Rscript cluster.R in.bed out.bed', call.=FALSE) }

df           <- read.table(args[1], header=FALSE)
colnames(df) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
df$cluster   <- kmeans(df$score,2)$cluster
#df$score     <- rescale(df$score,to=c(0,1000))

m <-df$cluster[order(-df$score)][1]
df$cluster[df$cluster!=m] <- 0
df$cluster[df$cluster==m] <- 1

df$score   <- df$cluster * df$score
df$cluster <- NULL

write.table(df[df$score!=0.0,], file=args[2], row.names=FALSE, quote=FALSE, col.names = FALSE)

# chr start end CNAG height strand (.bed columns)
