#install.packages("metap")
library(metap)

args <- commandArgs(trailingOnly = TRUE)
#args <- c('~/Projects/fixpeaks/pval-peak-by-gene/L1_OD1-pval.bed','~/Projects/fixpeaks/pval-peak-by-gene/L2_OD1-pval.bed','~/Projects/fixpeaks/pval-peak-by-gene/L_OD1-pval.bed')

rep1 <- read.table(args[1],header=FALSE,sep=" ")
rep2 <- read.table(args[2],header=FALSE,sep=" ")

names(rep1) <- c('chr','start','end','CNAG','pval','strand')
names(rep2) <- c('chr','start','end','CNAG','pval','strand')
rep1$pval   <- as.numeric(as.character(rep1$pval))
rep2$pval   <- as.numeric(as.character(rep2$pval))
rep1$pval   <- 10^(-rep1$pval)
rep2$pval   <- 10^(-rep2$pval)
rep1$pval[is.na(rep1$pval)] <- 1
rep2$pval[is.na(rep2$pval)] <- 1

rep1$corrected <- NA

n <- nrow(rep1)
for (i in 1:n){ rep1$corrected[i] <- sumlog(c(rep1$pval[i],rep2$pval[i]))$p }
rep1$corrected[rep1$corrected==1] <- NA
rep1$corrected <- -log10(rep1$corrected)
rep1$corrected[is.na(rep1$corrected)] <- '-'
rep1 <- rep1[,c(1,2,3,4,7,6)]

write.table(rep1,args[3],quote=FALSE,row.names=FALSE,col.names=FALSE,sep=" ")
