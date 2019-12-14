setwd('~/Projects/fixpeaks/edgeR-outputs/')

l <- read.table('L-edgeR-peaks.tsv',header=TRUE)
l <- na.omit(l)
View(l)
nrow(l)


ql <- read.table('qL-edgeR-peaks.tsv',header=TRUE)
ql <- na.omit(ql)
View(ql)
nrow(ql)

f <- 'N-edgeR-peaks.tsv'
n <- read.table(f,header=TRUE)
n <- na.omit(n)
View(n)
nrow(n)

f <- 'qN-edgeR-peaks.tsv'
qn <- read.table(f,header=TRUE)
qn <- na.omit(qn)
View(qn)
nrow(qn)

f <- 'C-edgeR-peaks.tsv'
c <- read.table(f,header=TRUE)
c <- na.omit(c)
View(c)
nrow(c)

f <- 'qC-edgeR-peaks.tsv'
qc <- read.table(f,header=TRUE)
qc <- na.omit(qc)
View(qc)
nrow(qc)

# ~~~~~~~~~~~~~~~~ differential ~~~~~~~~~~~~~~~~~
ld <- read.table('L-edgeR-diff.tsv',header=TRUE)
ld <- na.omit(ld)
nrow(ld)
View(ld)

cd <- read.table('C-edgeR-diff.tsv',header=TRUE)
cd <- na.omit(cd)
nrow(cd)
View(cd)

nd <- read.table('N-edgeR-diff.tsv',header=TRUE)
nd <- na.omit(nd)
nrow(nd)
View(nd)

# ~~~~~~~~ joint ~~~~~~~~
lld <- inner_join(ld,l,by='gene')
nrow(lld)
View(lld)
ccd <- inner_join(cd,c,by='gene')
nrow(ccd)

nnd <- inner_join(nd,n,by='gene')
nrow(nnd)

ngenes <- 6993

expected <- nrow(l) * nrow(n) * nrow(c) / ngenes^2
lc <- inner_join(l,c,by='gene')
lcn <- inner_join(lc,n,by='gene')
observed <- nrow(lcn)
observed / expected # enrichment

expected <- nrow(l) * nrow(n) / ngenes
ln <- inner_join(l,n,by='gene')
observed <- nrow(ln)
observed / expected

expected <- nrow(l) * nrow(c) / ngenes
lc <- inner_join(l,c,by='gene')
observed <- nrow(lc)
observed / expected

expected <- nrow(n) * nrow(c) / ngenes
nc <- inner_join(n,c,by='gene')
observed <- nrow(nc)
observed / expected






# ~~~~~~~~~ my method ~~~~~~~~~~~~~~
get_sum <- function(f){
c <- read.table(f)
names(c) <- c('chr','start','end','gene','strand','score')
print(head(c))
print(tail(c))
print(nrow(c))
return(c)
}

c  <- get_sum('C_OD1_3_enriched.bed')
qc <- get_sum('qC_OD1_3_enriched.bed')
n  <- get_sum('N_OD1_3_enriched.bed')
qn <- get_sum('qN_OD1_3_enriched.bed')
l  <- get_sum('L_OD1_3_enriched.bed')
ql <- get_sum('qL_OD1_3_enriched.bed')

525
709
780
2211
328
324


# ~~~~ enrichment fold change math with my method ~~~~~

ngenes <- 6993

expected <- nrow(l) * nrow(n) * nrow(c) / ngenes^2
lc <- inner_join(l,c,by='gene')
lcn <- inner_join(lc,n,by='gene')
observed <- nrow(lcn)
observed / expected # enrichment







