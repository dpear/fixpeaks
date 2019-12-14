# creating a table that is a summary of all values
library(extraDistr)
install.packages("d3heatmap")
library(d3heatmap)

# EDGE R FILES ~~~~~~~~~~~~~~~~
setwd('~/Projects/fixpeaks/edgeR-outputs/')
files <- c( dir('~/Projects/fixpeaks/edgeR-outputs/',pattern='*diff-all.tsv'),
            dir('~/Projects/fixpeaks/edgeR-outputs/',pattern='*diff.tsv'),
            dir('~/Projects/fixpeaks/edgeR-outputs/',pattern='*peaks-all.tsv'),
            dir('~/Projects/fixpeaks/edgeR-outputs/',pattern='*peaks.tsv')
)

df.whole <- read.table(files[1],header=TRUE)
df.whole <- na.exclude(df.whole)
colname   <- strsplit(files[1],'.tsv')[[1]][1]
names(df.whole) <- c('gene',
                     paste('logFC.',colname,sep=""),
                     paste('logCPM.',colname,sep=""),
                     paste('PValue.',colname,sep=""),
                     paste('FDR.',colname,sep=""))

for (f in files[2:length(files)]){
  df        <- read.table(f,header=TRUE)
  df        <- na.exclude(df)
  colname   <- strsplit(f,'.tsv')[[1]][1]
  names(df) <- c('gene',
                  paste('logFC.',colname,sep=""),
                  paste('logCPM.',colname,sep=""),
                  paste('PValue.',colname,sep=""),
                  paste('FDR.',colname,sep=""))
  df.whole <- full_join(df.whole,df,by='gene')
}


# ENRICHED SUMS FILES ~~~~~~~~~~~~~~~~~~~~~
setwd('~/Projects/fixpeaks/sums3-enriched/')
files2 <-  dir('~/Projects/fixpeaks/sums3-enriched/',pattern='[CNL][0-9]_OD1')
df.bed <- read.table(files2[1])
df.bed[,c(1,2,3,6)] <- NULL
colname   <- strsplit(files2[1],'.bed')[[1]][1]
names(df.bed) <- c('gene',paste('score.',colname,sep=''))

for (f in files2[2:length(files2)]){
  df <- read.table(f)
  df[,c(1,2,3,6)] <- NULL
  colname   <- strsplit(f,'.bed')[[1]][1]
  names(df) <- c('gene',paste('score.',colname,sep=''))
  df.bed <- full_join(df.bed,df,by='gene')
}

# combined replicates
files3 <- dir('~/Projects/fixpeaks/sums3-enriched/',pattern='[CNL]_OD1')
df.bed2 <- read.table(files3[1])
df.bed2[,c(1,2,3,5)] <- NULL
colname   <- strsplit(files3[1],'.bed')[[1]][1]
names(df.bed2) <- c('gene',paste('score.',colname,sep=''))

for (f in files3[2:length(files3)]){
  df <- read.table(f)
  df[,c(1,2,3,5)] <- NULL
  colname   <- strsplit(f,'.bed')[[1]][1]
  names(df) <- c('gene',paste('score.',colname,sep=''))
  df.bed2 <- full_join(df.bed2,df,by='gene')
}


# SUMS FILES ~~~~~~~~~~~~~~~~~~~~
setwd('~/Projects/fixpeaks/sums-peak-by-gene/')
files4 <- dir('~/Projects/fixpeaks/sums-peak-by-gene/',pattern='OD1')
df.bed3 <- read.table(files4[1])
df.bed3[,c(1,2,3,6)] <- NULL
colname   <- strsplit(files4[1],'.bed')[[1]][1]
names(df.bed3) <- c('gene',paste('score.',colname,sep=''))

for (f in files4[2:length(files4)]){
  df <- read.table(f)
  df[,c(1,2,3,6)] <- NULL
  colname   <- strsplit(f,'.bed')[[1]][1]
  names(df) <- c('gene',paste('score.',colname,sep=''))
  df.bed3 <- full_join(df.bed3,df,by='gene')
}

df.whole <- full_join(df.whole,df.bed,by='gene')
df.whole <- full_join(df.whole,df.bed3,by='gene')
df.whole <- full_join(df.whole,df.bed2,by='gene')


# RNA SEQ ~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('~/Projects/fixpeaks/')
df <- read.table('RNA_data_condensed.csv',sep=',',header = TRUE)
names(df)[1] <- 'gene'

df.whole <- full_join(df.whole,df,by='gene')


write.table(df.whole,'~/Projects/fixpeaks/whole-table.tsv',quote=FALSE,col.names=TRUE,row.names=FALSE)


# Diana's table

setwd('~/Projects/fixpeaks/')

C.file  <- c('sums-peak-by-gene/C1_OD1_sums.bed','sums-peak-by-gene/C2_OD1_sums.bed','sums3-enriched/C_OD1_3_enriched.bed')
#lC.file <- c('sums-peak-by-gene-14/lC1_sums.bed','sums-peak-by-gene-14/lC2_sums.bed','sums3-enriched-14/lC_3_enriched.bed')
nC.file <- c('sums-peak-by-gene-14/nC1_sums.bed','sums-peak-by-gene-14/nC2_sums.bed','sums3-enriched-14/nC_3_enriched.bed')
qC.file <- c('sums-peak-by-gene/qC1_OD1_sums.bed','sums-peak-by-gene/qC2_OD1_sums.bed','sums3-enriched/qC_OD1_3_enriched.bed')
C14.file <- c('sums-peak-by-gene-14/C1_sums.bed','sums-peak-by-gene-14/C2_sums.bed','sums3-enriched-14/C_3_enriched.bed')
  
cL.file <- c('sums-peak-by-gene-14/cL1-b_sums.bed','sums-peak-by-gene-14/cL2-b_sums.bed','sums3-enriched-14/cL-b_3_enriched.bed')
L.file  <- c('sums-peak-by-gene/L1_OD1_sums.bed','sums-peak-by-gene/L2_OD1_sums.bed','sums3-enriched/L_OD1_3_enriched.bed')
nL.file <- c('sums-peak-by-gene-14/nL1_sums.bed','sums-peak-by-gene-14/nL2_sums.bed','sums3-enriched-14/nL_3_enriched.bed') 
qL.file <- c('sums-peak-by-gene/qL1_OD1_sums.bed','sums-peak-by-gene/qL2_OD1_sums.bed','sums3-enriched/qL_OD1_3_enriched.bed')  
L14.file <- c('sums-peak-by-gene-14/L1_sums.bed','sums-peak-by-gene-14/L2_sums.bed','sums3-enriched-14/L_3_enriched.bed')

cN.file <- c('sums-peak-by-gene-14/cN1_sums.bed','sums-peak-by-gene-14/cN2_sums.bed','sums3-enriched-14/cN_3_enriched.bed')
#lN.file <- c('sums-peak-by-gene-14/lN1_sums.bed','sums-peak-by-gene-14/lN2_sums.bed','sums3-enriched-14/lN_3_enriched.bed')
N.file  <- c('sums-peak-by-gene-14/N1_sums.bed' ,'sums-peak-by-gene-14/N2_sums.bed' ,'sums3-enriched-14/N_3_enriched.bed')
qN.file <- c('sums-peak-by-gene-14/qN1_sums.bed','sums-peak-by-gene-14/qN2_sums.bed','sums3-enriched-14/qN_3_enriched.bed')

all.files <- cbind(C.file,nC.file,qC.file,C14.file,cL.file,L.file,nL.file,qL.file,L14.file,cN.file,N.file,qN.file)

samples <- c('C','nC','qC','C14','cL','L','nL','qL','L14','cN','N','qN')

es <- 1:length(samples)

full <- data.frame(gene='cnag')
for (e in es){

  sample <- samples[e]

  r1 <- read.table(all.files[,e][1])
  r2 <- read.table(all.files[,e][2])
  bound <- read.table(all.files[,e][3])
  
  r1[,c(1,2,3,6)] <- NULL
  r2[,c(1,2,3,6)] <- NULL
  bound[,c(1,2,3,5)] <- NULL
  
  #colnames <- c('gene',paste(sample,'.score',sep=''))
  colnames <- c('gene','score')
  names(r1) <- colnames
  names(r2) <- colnames
  names(bound) <- colnames
  
  reps <- full_join(r1,r2,     by='gene')
  all  <- full_join(reps,bound,by='gene')
  all$bound <- as.numeric(!is.na(all$score))
  all$score <- NULL
  all$score <- apply(all[,2:3],1,max)
  all[,2:3] <- NULL
  colnames <- c('gene',paste(sample,'.bound',sep=''),paste(sample,'.score',sep=''))
  names(all) <- colnames
  
  full <- full_join(full,all,by='gene')
}
full <- full[2:nrow(full),]

wt.cols <- c('C14.score','C.score','L14.score','L14.score','L.score','N.score','N.score')
mt.cols <- c('nC.score','qC.score','cL.score','nL.score','qL.score','cN.score','qN.score')
wt.bound <- c('C14.bound','C.bound','L14.bound','L14.bound','L.bound','N.bound','N.bound')
mt.bound <- c('nC.bound','qC.bound','cL.bound','nL.bound','qL.bound','cN.bound','qN.bound')
new.cols <- c('C14.nC.fold','C.qC.fold','L14.cL.fold','L14.nL.fold','L.qL.fold','N.cN.fold','N.qN.fold')
diff.cols <- c('C14.nC.diff','C.qC.diff','L14.cL.diff','L14.nL.diff','L.qL.diff','N.cN.diff','N.qN.diff')
cutoff <- 1.5
full[,new.cols] <- full[,mt.cols]/full[,wt.cols]
full[,diff.cols] <- (full[,new.cols] > cutoff | full[,new.cols] < 1/cutoff) & (full[,wt.bound] | full[,mt.bound])
full[,diff.cols] <- full[,diff.cols] * 1
full[,new.cols] <- log2(full[,new.cols])




head(full)

names(full)

# RNA DATA

rna.files <- c('/data/diana/BRHM0304_analysis/WT_v_C1_YNB_DE.csv',  
               '/data/diana/BRHM0304_analysis/WT_v_L1_YNB_DE.csv',  
               '/data/diana/BRHM0304_analysis/WT_v_N1_YNB_DE.csv',
               '/data/diana/BRHM0304_analysis/WT_v_Q1_YNB_DE.csv',
               '/data/diana/BRHM0304_analysis/WT_v_C1_YNB_DE_sig_changed.csv',
               '/data/diana/BRHM0304_analysis/WT_v_L1_YNB_DE_sig_changed.csv',
               '/data/diana/BRHM0304_analysis/WT_v_N1_YNB_DE_sig_changed.csv',
               '/data/diana/BRHM0304_analysis/WT_v_Q1_YNB_DE_sig_changed.csv',
               '/data/diana/BRHM0304_analysis/WT_v_C1_YNB_DE_increased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_L1_YNB_DE_increased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_N1_YNB_DE_increased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_Q1_YNB_DE_increased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_C1_YNB_DE_decreased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_L1_YNB_DE_decreased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_N1_YNB_DE_decreased.csv',
               '/data/diana/BRHM0304_analysis/WT_v_Q1_YNB_DE_decreased.csv')

samples <- c('WT_v_C1_all',  
             'WT_v_L1_all',
             'WT_v_N1_all',
             'WT_v_Q1_all',
             'WT_v_C1_sig_changed',
             'WT_v_L1_sig_changed',
             'WT_v_N1_sig_changed',
             'WT_v_Q1_sig_changed',
             'WT_v_C1_increased',
             'WT_v_L1_increased',
             'WT_v_N1_increased',
             'WT_v_Q1_increased',
               'WT_v_C1_decreased',
               'WT_v_L1_decreased',
               'WT_v_N1_decreased',
               'WT_v_Q1_decreased')

readin <- function(df,file,sample){
  rna.table <- read.csv(file)
  rna.table[,c(2,4:7)] <- NULL
  names(rna.table) <- c('gene',sample)
  df <- full_join(df,rna.table,by='gene')
  return(df)
}

for (i in 1:length(samples)){ full <- readin(full,rna.files[i],samples[i])}



rna.cols <- names(full)[44:55]
full[,rna.cols] <- !is.na(full[,rna.cols])
full[,rna.cols] <- full[,rna.cols]*1

# remove cnag_1****

full <- full[!grepl('CNAG_1',full$gene),]



head(full)

# Just RNA table

just_rna <- semi_join(full,f,by='gene')
just_rna <- just_rna[,names(just_rna)[c(1,2,3,6,7,12,13,16,17,22,23,24,25,27,30,32,34,37,39,40,41,42,43)]]
View(just_rna)
write.table(just_rna,file='~/Projects/fixpeaks/just_rna_data.csv',quote=FALSE,row.names=FALSE,sep=',')








# HEATMAP

# genes that are differentially bound and must be bound by C L or N
any.diff <-    (full$C.qC.diff==1 | full$L.qL.diff==1 | full$N.qN.diff==1)  & 
               (full$C.bound==1 | full$L.bound==1 | full$N.bound==1 | full$qC.bound | full$qL.bound | full$qN.bound )
C.diff <- full$C.qC.diff==1
L.diff <- full$L.qL.diff==1
N.diff <- full$N.qN.diff==1
rna.diff <- full$WT_v_Q1_sig_changed==1

names(full)
hm <- full[ any.diff
            ,
          c('gene','C.qC.fold','L.qL.fold','N.qN.fold','WT_v_Q1_all',
                   'C.qC.diff','L.qL.diff','N.qN.diff','WT_v_Q1_sig_changed')]
hm <- hm[complete.cases(hm[ , 1]),]
#hm[hm==NA] <- 0

values <- as.matrix(hm[,2:5]*hm[,6:9])
values <- as.data.frame(values)
values[is.na(values)] <- 0




d3heatmap(as.matrix(values),
          labRow=hm[,1],
          labCol=c('Cqs2','Liv3','Nrg1','WT vs. Q1'),
          cexCol = .8,
          na.color='white',
          #Rowv = FALSE,
          #Colv = FALSE,
          scale = "none",
          dendrogram = 'none',
          #col=colorRampPalette(c("purple", "blue", "green","yellow")),
          #trace = 'none'
          #col='virids'
)

map <- heatmap.2(as.matrix(values),
        labRow=hm[,1],
        labCol=c('Cqs2','Liv3','Nrg1','WT vs. Q1'),
        cexCol = .8,
        na.color='white',
        # Rowv = FALSE,
        # Colv = FALSE,
        scale = "none",
        #dendrogram = 'none',
        col=colorRampPalette(c("purple", "blue", "green","yellow")),
        #trace = 'none'
        #col='virids'
        )



setwd('../edgeR-outputs/')
diff.files <- dir('.',pattern='diff-all')
all.files <- dir('.',pattern='diff.tsv')

for (i in 1:length(all.files)){
  
  sample <- strsplit(diff.files[i],'-')[[1]][1]
  
  diff <- read.table(diff.files[i],header=TRUE)
  all <- read.table(all.files[i],header=TRUE)
  
  diff[,3:4] <- NULL
  all[,3:5] <- NULL
  
  colname <- c('gene',paste(sample,'.edgeR.logFC',sep=''),paste(sample,'.edgeR.FDR',sep=''))
  colnames(diff) <- colname
  colnames(all) <- colname[1:2]
  colnames(all) <- c('gene',paste(sample,'.edgeR.diff',sep=''))
  all <- full_join(all,diff,by='gene')
  
  
  all[,2] <- as.numeric(!is.na(all[,2]))
  full <- full_join(full,all,by='gene')
  colnames(all)[2] <- paste(sample,'.edgeR.bound',sep='')
                     
}

names(full)
a <- na.omit(full[,c('gene','C.logFC')][full$C.bound==1 & full$C.logFC > .4 ,])
b <- na.omit(full[,c('gene','C.edgeR.logFC')][full$C.bound==1 & full$C.edgeR.logFC < -.4, ])
c <- inner_join(a,b,by='gene')
c$C.logFC <- -c$C.logFC
head(c)

# Statistical analysis
N <- nrow(full)

f1 <- 'N.bound'
f2 <- 'cN.bound'

hyper.test <- function(f1,f2){
  f1.table <- na.omit(full[,c('gene',f1)][full[,f1]==1,])
  f2.table <- na.omit(full[,c('gene',f2)][full[,f2]==1,])
  f.table <- inner_join(f1.table, f2.table, by='gene')
  
  q <- nrow(f.table)
  m <- nrow(f1.table)
  n <- N-m
  k <- nrow(f2.table)
  pval <- 1-phyper(q,m,n,k)
  return(c(m,k,q,pval))
}


N.cN <- hyper.test('N.bound','cN.bound')
L.cL <- hyper.test('L14.bound','cL.bound')
L.nL <- hyper.test('L14.bound','nL.bound')
cL.nL <- hyper.test('cL.bound','nL.bound')
C.nC <- hyper.test('C14.bound','nC.bound')

N.qN_N.cN <- hyper.test('N.qN.diff','N.cN.diff')
#N.qN_N.lN <- hyper.test()
#N.cN_N.lN <- hyper.test()
#C.lC_C.nC <- hyper.test()
#C.lC_C.qC <- hyper.test()
C.nC_C.qC <- hyper.test('C14.nC.diff','C.qC.diff')
L.cL_L.nL <- hyper.test('L14.cL.diff','L14.nL.diff')
L.cL_L.qL <- hyper.test('L14.cL.diff','L.qL.diff')
L.nL_L.qL <- hyper.test('L14.nL.diff','L.qL.diff')



N.qN <- hyper.test('N.bound','qN.bound')
C.qC <- hyper.test('C.bound','qC.bound')
L.qL <- hyper.test('L.bound','qL.bound')


# ChIP/RNA comparisons

C.WT_v_C1_sig_changed <- hyper.test('C.bound','WT_v_C1_sig_changed')
L.WT_v_L1_sig_changed <- hyper.test('L.bound','WT_v_L1_sig_changed')
N.WT_v_N1_sig_changed <- hyper.test('N.bound','WT_v_N1_sig_changed')
C_qC.WT_v_Q1_sig_changed <- hyper.test('C.qC.diff','WT_v_Q1_sig_changed')
L_qL.WT_v_Q1_sig_changed <- hyper.test('L.qL.diff','WT_v_Q1_sig_changed')
N_qN.WT_v_Q1_sig_changed <- hyper.test('N.qN.diff','WT_v_Q1_sig_changed')

C.WT_v_C1_increased <- hyper.test('C.bound','WT_v_C1_increased')
L.WT_v_L1_increased <- hyper.test('L.bound','WT_v_L1_increased')
N.WT_v_N1_increased <- hyper.test('N.bound','WT_v_N1_increased')
C_qC.WT_v_Q1_increased <- hyper.test('C.qC.diff','WT_v_Q1_increased')
L_qL.WT_v_Q1_increased <- hyper.test('L.qL.diff','WT_v_Q1_increased')
N_qN.WT_v_Q1_increased <- hyper.test('N.qN.diff','WT_v_Q1_increased')

C.WT_v_C1_decreased <- hyper.test('C.bound','WT_v_C1_decreased')
L.WT_v_L1_decreased <- hyper.test('L.bound','WT_v_L1_decreased')
N.WT_v_N1_decreased <- hyper.test('N.bound','WT_v_N1_decreased')
C_qC.WT_v_Q1_decreased <- hyper.test('C.qC.diff','WT_v_Q1_decreased')
L_qL.WT_v_Q1_decreased <- hyper.test('L.qL.diff','WT_v_Q1_decreased')
N_qN.WT_v_Q1_decreased <- hyper.test('N.qN.diff','WT_v_Q1_decreased')

# checks

C.C <- hyper.test('C.bound','C14.bound')






















# RNA DATA
setwd('~/Projects/fixpeaks/')
f <- read.table('RNA_data_condensed.csv',header = TRUE,sep=',')
names(f)[1] <- 'gene'
head(f)
full <- full_join(full,f,by='gene')


head(full)

write.table(full,file='~/Projects/fixpeaks/whole-data-table-new.csv',quote=FALSE,row.names=FALSE,sep=',')



nrow(full[  full$C.bound==1 & full$N.bound==1 & full$L.bound==1,]) #310 for all 3

# binomial test
run.binom.test <- function(tf1,tf2,conservative=FALSE){
  n=6993
  if(conservative){n=nrow(full[full$C.bound==1 | full$L.bound==1 | full$N.bound==1,])}
  if (tf1=='c' & tf2=='l'){
    k <- nrow(full[  full$C.bound==1 & full$L.bound==1,])
    p <- nrow(full[full$C.bound==1,])*nrow(full[full$L.bound==1,])/6993^2
    b <- binom.test(k,n,p)
  }
  
  if (tf1=='c' & tf2=='n'){
    k <- nrow(full[  full$C.bound==1 & full$N.bound==1,])
    p <- nrow(full[full$C.bound==1,])*nrow(full[full$N.bound==1,])/6993^2
    b <- binom.test(k,n,p)
  }
  
  if (tf1=='l' & tf2=='n'){
    k <- nrow(full[  full$L.bound==1 & full$N.bound==1,])
    p <- nrow(full[full$L.bound==1,])*nrow(full[full$N.bound==1,])/6993^2
    b <- binom.test(k,n,p)
  }
  
  if(tf1=='all' & tf2=='all'){
    k <- nrow(full[  full$L.bound==1 & full$N.bound==1 & full$C.bound==1,])
    p <- nrow(full[full$L.bound==1,])*nrow(full[full$N.bound==1,])*nrow(full[full$C.bound==1,])/6993^3
    b <- binom.test(k,n,p)
  }
  
  return(b)
  
}

b_cl <- run.binom.test('c','l')
b_cn <- run.binom.test('c','n')
b_ln <- run.binom.test('l','n')
b_all<- run.binom.test('all','all') 



run.hypergeometric.test <- function(tf1,tf2,conservative=FALSE){
  
  N <- 6993
  if(conservative){ N <- length( na.omit(full$C.bound[full$C.bound==1 | full$L.bound==1 | full$N.bound==1]) )  }
  
  if (tf1=='c' & tf2=='l'){
    m <- length( na.omit(full$C.bound[full$C.bound==1]) )# all tf1 genes
    n <- N-m
    k <- length( na.omit(full$C.bound[full$L.bound==1]) )# all tf2 genes
    x <- length( na.omit(full$C.bound[full$C.bound==1 & full$L.bound==1]) )# tf1 and tf2 genes
    h <- 1-phyper(x,m,n,k)
  }
  
  if (tf1=='c' & tf2=='n'){
    m <- length( na.omit(full$C.bound[full$C.bound==1]) )# all tf1 genes
    n <- N-m
    k <- length( na.omit(full$C.bound[full$N.bound==1]) )# all tf2 genes
    x <- length( na.omit(full$C.bound[full$C.bound==1 & full$N.bound==1]) )# tf1 and tf2 genes
    h <- 1-phyper(x,m,n,k)
  }
  
  if (tf1=='l' & tf2=='n'){
    m <- length( na.omit(full$C.bound[full$L.bound==1]) )# all tf1 genes
    n <- N-m
    k <- length( na.omit(full$C.bound[full$N.bound==1]) )# all tf2 genes
    x <- length( na.omit(full$C.bound[full$L.bound==1 & full$N.bound==1]) )# tf1 and tf2 genes
    h <- 1-phyper(x,m,n,k)
  }
  
  if(tf1=='all' & tf2=='all'){

    o   <- length( na.omit(full$C.bound[full$C.bound==0 & full$L.bound==0 & full$N.bound==0]) )
    c   <- length( na.omit(full$C.bound[full$C.bound==1 & full$L.bound==0 & full$N.bound==0]) )
    l   <- length( na.omit(full$C.bound[full$C.bound==0 & full$L.bound==1 & full$N.bound==0]) )
    n   <- length( na.omit(full$C.bound[full$C.bound==0 & full$L.bound==0 & full$N.bound==1]) )
    cl  <- length( na.omit(full$C.bound[full$C.bound==1 & full$L.bound==1 & full$N.bound==0]) )
    cn  <- length( na.omit(full$C.bound[full$C.bound==1 & full$L.bound==0 & full$N.bound==1]) )
    ln  <- length( na.omit(full$C.bound[full$C.bound==0 & full$L.bound==1 & full$N.bound==1]) )
    cln <- length( na.omit(full$C.bound[full$C.bound==1 & full$L.bound==1 & full$N.bound==1]) )
    
    ni <- c(c+cn,o+n,l+ln,cl+cln)
    xi <- c(cn,n,ln,cln)
    k <- sum(xi)
    
    # using a multinomial hypergeometric distribution: ?dmvhyper
    h <- dmvhyper(x=xi,n=ni,k=k)
  }
  
  return(h)
}

h_cl <- run.hypergeometric.test('c','l',conservative = TRUE)
h_cn <- run.hypergeometric.test('c','n')
h_ln <- run.hypergeometric.test('l','n')
h_all<- run.hypergeometric.test('all','all')


# venn diagram

c <- na.omit(full$gene[full$C.bound==1])
n <- na.omit(full$gene[full$N.bound==1])
l <- na.omit(full$gene[full$L.bound==1])
genes <- list(c,n,l)
venn.diagram(x = genes,
             category.names = c('c','n','l'),
             filename = 'venn-diagram.png',
             output = TRUE)




# hypergeometric test
g <- 75 ## Number of submitted genes
k <- 59 ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
m <- 611 ## Number of "marked" elements, i.e. genes associated to this biological process
N <- 13588 ## Total number of genes with some annotation in GOTERM_BP_FAT.  
n <- N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
x <- 19 ## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process


m <- nrow(full[full$C.bound==1 | full$N.bound==1,])
N <- 6993
n <- N-m
k <- nrow(full[full$C.bound==1 & full$N.bound==1,])
x <- 

p.value <-  phyper(q=x -1, m=m, n=n, k=k)



scatterlist <- list()

data <- full[full$C.bound==1 | full$qC.bound==1,][complete.cases(full[full$C.bound==1 | full$qC.bound==1,1]),]

p <- ggplot(xlab='score') + 
  geom_point(data=data, aes(x=C.score, y=qC.score, color=C.qC.diff)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(title='Nrg1', x='WT score',y='Mutant score') +
  scale_x_log10( limits=c(1000,80000)) + 
  scale_y_log10( limits=c(1000,80000)) +
  coord_fixed() +
  theme(axis.text.x = element_text(face="italic", color="black", 
                                   size=6, angle=45),
        axis.text.y = element_text(face="italic", color="black", 
                                   size=6, angle=45),
        panel.background = element_rect(fill = "white",
                                        colour = "grey",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey")) 
  
  
  

p
scatterlist[[1]] <- p

plot_scatters <- function(df,title,subtitle1,subtitle2){   
  sc <- ggplot(df, aes(x=wt_score,y=mut_score)) + 
    geom_point(alpha = 0.18, aes(x=wt_score, y=mut_score, color=color), show.legend = FALSE) +
    geom_abline(slope=1,intercept=0,col="grey") +
    scale_x_log10( limits=c(1000,80000)) + 
    scale_y_log10( limits=c(1000,80000)) +
    labs( title=title, subtitle=subtitle1, x=subtitle2,y="") +
    theme(axis.text.x = element_text(face="italic", color="black", 
                                     size=6, angle=45),
          axis.text.y = element_text(face="italic", color="black", 
                                     size=6, angle=45),
          panel.background = element_rect(fill = "white",
                                          colour = "grey",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey")) +
    scale_color_manual(values=c("darkblue","red")) +
    coord_fixed()  # force plot boundaries to remain square
  #theme( panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"))
  
  s <- sc #+ ggtitle(title,subtitle=sub)
  return(s)
}

length(na.omit(full$C.score[full$C.bound==1 & full$C.edgeR.diff==1]))
length(na.omit(full$C.score[full$N.bound==1 & full$N.edgeR.diff==1]))
length(na.omit(full$C.score[full$L.bound==1 & full$L.edgeR.diff==1]))

"""
p <- ggplot() + 
  geom_point(data=full, aes(x=eval(parse(text=dataset[1])), y=eval(parse(text=dataset[2])), color=eval(parse(text=dataset[3])))) +
  geom_abline(slope = 1, intercept = 0) + 
  
  p"""


grid.arrange(
  scatterlist[[1 ]],scatterlist[[2 ]],scatterlist[[3 ]],
  ncol=3,nrow=1,
  widths = c(1,1,1),
  top = 'WT vs. Mutant Scores (Using Sums3) Showing bound genes only',
  left ="Mutanat Score",
  bottom="Wild Type Score"
)


  
f# change diff (not diff-all) to 1 or 0

no.cpm <- !grepl('logCPM',colnames(df.whole))
no.pval <- !grepl('PValue',colnames(df.whole))
no.fc.diff <- !( grepl('logFC',colnames(df.whole)) & grepl('diff$',colnames(df.whole)))
no.fc.peaks <- !( grepl('logFC',colnames(df.whole)) & grepl('peaks$',colnames(df.whole)))

logfc.diff.all <- grepl('logFC',colnames(df.whole)) & grepl('diff-all',colnames(df.whole))
genes <- grepl('gene',colnames(df.whole))

condition <- logfc.diff.all | genes
colnames(df.whole[,condition])

heat.df <- df.whole[,condition]
colnames(heat.df) <- c('gene','C','L','N')
colnames(heat.df)


sum(is.infinite(as.matrix(heat.df[,2:ncol(heat.df)])))
#heat.df[is.na(heat.df)] <- 10
#heat.df <- heat.df[rowSums(heat.df[2:ncol(heat.df)])!=30,]
#heat.df <- heat.df[rowSums(heat.df[2:ncol(heat.df)])< -2, ]
heat.df <- heat.df[order(heat.df$C),]
#heat.df <- head(heat.df,1000)
heatmap.2(as.matrix(heat.df[2:ncol(heat.df)]),scale='none',dendogram='row',labRow = heat.df$gene, na.color = "white",Rowv = FALSE, Colv=FALSE)

# preset conditions
logFC.ind    <- grepl('logFC',colnames(df.whole))
FDR.ind      <- grepl('FDR',colnames(df.whole))
C            <- grepl('C-',colnames(df.whole))
peaks        <- grepl('peaks',colnames(df.whole))
notall       <- !grepl('all',colnames(df.whole))

these.columns <- df.whole[,c('PValue.C-edgeR-diff-all',"PValue.C-edgeR-diff")]


condition <- these.columns
View(head(df.whole[,condition]))

test <- head(df.whole[,condition])


is.na(test$`FDR.qC-edgeR-peaks`)
