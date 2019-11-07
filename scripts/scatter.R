# demo(plotmath) # math symbols in R
# fancy scientific didn't work for labels

files = list.files('~/Projects/fixpeaks/peak-by-gene-outputs',full.names = TRUE)

wt <- read.table(files[1],header=FALSE,sep="\t")
mut <- read.table(files[12],header=FALSE,sep="\t")
wt <- wt[order(wt$V4),]
mut <- mut[order(mut$V4),]

df <- data.frame(cnag=wt$V4,wt=wt$V5,mut=mut$V5)
plt <- ggplot(df, aes(x=wt-mut)) +
  geom_density(alpha=0.7)


enrich_files <- list.files("~/Projects/fixpeaks/sums-enriched/",full.names = TRUE)
mut_files <- list.files("~/Projects/fixpeaks/sums-peak-by-gene/",full.names = TRUE)

plotlist <- list()
scatterlist <- list()
all <- data.frame(a=NA,b=NA,c=NA,CNAG=NA,wt_score=NA,g=NA,d=NA,e=NA,f=NA,mut_score=NA,h=NA,combined=NA,color=NA,name=NA)
#percent_fold <- list()
n=length(enrich_files)/2
diff <- data.frame(CNAG=NA, wt_score=NA, mut_score=NA, combined=NA, color=NA, name=NA)

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

plot_scatters(joined,title,subtitle,subtitle2)

for (i in 1:n){
  enriched_wt        <- read.table(enrich_files[i],header=FALSE,sep=" ")
  names(enriched_wt) <- c('a','b','c','CNAG','wt_score','g')

  all_mut         <- read.table(mut_files[i+n])
  names(all_mut)  <- c('d','e','f','CNAG','mut_score','h')

  joined          <- left_join(enriched_wt,all_mut)
  joined$combined <- joined$wt_score/joined$mut_score
  joined$color    <- 'darkblue'
  joined$color[ joined$combined>2 | joined$combined<.5 ] <- 'red'
  
  title        <- strsplit(strsplit(enrich_files[i],"//")[[1]][2], '_e')[[1]][1]
  total        <- length(joined$combined)
  
  num          <- sum(joined$combined>1.5)
  percent_fold <- round(num/total,4)
  subtitle     <- paste(paste(num,"/",total,sep=""),"genes: FC > 2 (",percent_fold,")")
  
  num2          <- sum(joined$combined<.5)
  percent_fold2 <- round(num2/total,4)
  subtitle2     <- paste(paste(num2,"/",total,sep=""),"genes: FC < .5 (",percent_fold2,")")
  
  plotlist[[i]]    <- plot_multi_histogram(joined,title)
  scatterlist[[i]] <- plot_scatters(joined,title,subtitle,subtitle2)
  
  different <- joined[joined$combined>1.5,]
  different[,c(1,2,3,6,7,8,9,11)] <- NULL
  different$name <- title
  diff <- rbind(diff,different)
  joined$name <- title
  all <- rbind(all,joined)
}

grid.arrange(
  scatterlist[[1 ]],scatterlist[[8 ]],scatterlist[[5 ]],
  scatterlist[[3 ]],scatterlist[[10]],scatterlist[[7 ]],
  scatterlist[[2 ]],scatterlist[[9 ]],scatterlist[[6 ]],
  scatterlist[[4 ]],scatterlist[[11]],ggplot(),
  ncol=3,nrow=4,
  widths = c(1,1,1),
  top = 'WT vs. Mutant Scores (Using Sums)',
  left ="Mutanat Score",
  bottom="Wild Type Score"
)

diff <- na.omit(diff)

joined <- inner_join(diff[diff$name=='C1_OD1',],diff[diff$name=='C2_OD1',],by='CNAG')
joined
fout <- '~/Projects/fixpeaks/all_fc_5.txt'
write.table(joined,fout,quote=FALSE,sep=" ",row.names = FALSE, col.names = TRUE)


plot_multi_histogram <- function(df ,title) {
  plt <- ggplot(df, aes(x=combined))+
    geom_density(alpha=0.7, fill="black") + coord_cartesian( xlim = c(0,10)) + 
    geom_vline(aes(xintercept=1), linetype="dashed",col="red") +
    labs(x="", y = "")
  p <- plt + ggtitle(title) 
  return(p)
}


grid.arrange(
  plotlist[[1 ]],plotlist[[8 ]],plotlist[[5 ]],
  plotlist[[3 ]],plotlist[[10]],plotlist[[7 ]],
  plotlist[[2 ]],plotlist[[9 ]],plotlist[[6 ]],
  plotlist[[4 ]],plotlist[[11]],ggplot(),
  ncol=3,nrow=4,
  widths = c(.5,.5,.5),
  top = 'Ratio of WT vs. Mutant',
  left ="Number of Genes",
  bottom="Fold Change"
)


ggplot() + labs(title = "Title size", subtitle = expression(atop(paste(.(subtitle)), subtitle2)))
assay <- subtitle
xlab <- bquote(.(assay))

plot(0, xlab = subtitle)
