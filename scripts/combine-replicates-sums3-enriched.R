setwd('~/Projects/fixpeaks/sums3-enriched-14/')


args <- commandArgs(trailingOnly = TRUE)
args <- c('C1_3_enriched.bed','C2_3_enriched.bed','C_3_enriched.bed')

rep1 <- read.table(args[1])
rep2 <- read.table(args[2])

names(rep1) <- c('chr','start','end','gene','score','strand')
names(rep2) <- c('chr','start','end','gene','score','strand')

both <- rbind(rep1,rep2)
both <- aggregate(both$score,by=list(chr=both$chr,
                               start=both$start,
                               end=both$end,
                               gene=both$gene,
                               strand=both$strand),max)

both <- semi_join(both,rep1,by='gene')
both <- semi_join(both,rep2,by='gene')
both <- both[order(-both$x),]


write.table(both,args[3],row.names = FALSE,col.names = FALSE,quote=FALSE)
