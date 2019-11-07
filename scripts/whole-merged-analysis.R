install.packages("OutlierDetection")



mergefile="~/Projects/fixpeaks/merged"
m = read.table(mergefile,sep=" ",header=TRUE)
View(m)

folder = "~/Projects/fixpeaks/peak-by-gene-outputs/"
ranks = list.files(folder, full.names = TRUE)

read_all <- function(f){ 
  tab = read.table(f,sep ="\t" ,header=FALSE)
  colnames(tab) <- c('chr','start','end','cnag','height','strand')
  return(tab)
}


all <- lapply(ranks,read_all)
names(all) <- list.files(folder, full.names = FALSE)

max(abs(m$C1_OD1_height))# - m$C2_OD1_height))

plot(all$C1_OD5_final.bed$,all$C2_OD5_final.bed$height)

par(mfrow=c(1,1))
par(mfrow=c(4,3))
plot(m$C1_OD1_height,m$C2_OD1_height,main="C_OD1")
plot(m$N1_OD1_height,m$N2_OD1_height,main="N_OD1")
plot(m$L1_OD1_height,m$L2_OD1_height,main="L_OD1")
plot(m$C1_OD5_height,m$C2_OD5_height,main="C_OD5")
plot(m$N1_OD5_height,m$N2_OD5_height,main="N_OD5")
plot(1,1,main="L_OD5")
plot(m$qC1_OD1_height,m$qC2_OD1_height,main="qC_OD1")
plot(m$qN1_OD1_height,m$qN2_OD1_height,main="qN_OD1")
plot(m$qL1_OD1_height,m$qL2_OD1_height,main="qL_OD1")
plot(m$qC1_OD5_height,m$qC2_OD5_height,main="qC_OD5")
plot(1,1,main="qN_OD5")
plot(1,1,main="qL_OD5")

head(m$C1_OD1_height[order(m$C1_OD1_height)])

var(m$C1_OD1_height,m$C2_OD1_height,main="C_OD1")


samp=m$C2_OD5_height
plot(samp[order(-samp)])

# sorted height plots for all samples
l=length(all)
for (i in 1:l){
  plot(all[[i]]$height[order(-all[[i]]$height)])
}


cidx=seq(7,47,by=2)

for (i in cidx){
  m[,i-1]=kmeans(m[,i],2)$cluster
}

for (i in cidx){
  plot(m[,i][order(-m[,i])],col=m[,i-1][order(-m[,i])])
}

colnames(m)[7]
head(m[,7])
paste("cluster",colnames(m)[7])

fit <- kmeans(s,2)


