library(ggplot2)
library(gridExtra)
setwd("~/Projects/fixpeaks/enriched/")


#args <- c('~/Projects/fixpeaks/enriched/N2_OD1_enriched.bed','~/Projects/fixpeaks/enriched/qN2_OD1_enriched.bed')
files <- list.files('~/Projects/fixpeaks/sums-enriched/',full.names = TRUE)
n <- length(files)
df <- data.frame(chr = NA, start=NA, end=NA, CNAG=NA, score=NA, strand=NA, condition=NA, replicate=NA, sample=NA, od=NA, tf=NA)
for (i in 1:n){
  sample          <- strsplit(strsplit(files[i],'_e')[[1]][1],'//')[[1]][2]
  genes           <- read.table(files[i],header=FALSE,sep=" ")
  names(genes)    <- c('chr','start','end','CNAG','score','strand')
  
  if     ( substr(sample,1,1)=='q'){ 
           condition <- 'Mutant'
  } else { condition <- 'WildType' }
  
  if     ( condition=='Mutant'){
           od <- paste('OD',substr(sample,7,7),sep='')
  } else { od <- paste('OD',substr(sample,6,6),sep='')}
  
  if     ( condition=='Mutant'){
           replicate <- substr(sample,3,3)
  } else { replicate <- substr(sample,2,2)}
   
  if     ( condition=='Mutant'){
           letter <- substr(sample,2,2)
  } else { letter <- substr(sample,1,1)}
  if (letter=='C'){ tf <- 'Cqs2'}
  if (letter=='N'){ tf <- 'Nrg1'}
  if (letter=='L'){ tf <- 'Liv3'}
  
  
  genes$condition <- condition
  genes$replicate <- replicate
  genes$sample    <- sample
  genes$od        <- od
  genes$tf        <- tf
  df              <- rbind(df,genes)
}

df <- na.omit(df)

df$in_both_replicates <- '_'
df$in_both_replicates <- aggregate()

#testing purposes
littledf <- df[df$CNAG=="CNAG_08008", ]

