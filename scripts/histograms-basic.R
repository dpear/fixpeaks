library(ggplot2)
library(gridExtra)
setwd("~/Projects/fixpeaks/enriched/")

plot_multi_histogram <- function(df, feature, label_column,title) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_density(alpha=0.7) + #coord_cartesian( ylim = c(0,.009)) + 
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), linetype="dashed", size=.5) +
    labs(x=feature, y = "")
  p <- plt + guides(fill=guide_legend(title=label_column)) + ggtitle(title,) 
  return(p)
}

#args <- c('~/Projects/fixpeaks/enriched/N2_OD1_enriched.bed','~/Projects/fixpeaks/enriched/qN2_OD1_enriched.bed')
files <- list.files('~/Projects/fixpeaks/sums-enriched/',full.names = FALSE)
n <- length(files)/2
plotlist <- list()
for (i in 1:n){
  genes <- read.table(files[i],header=FALSE,sep=" ")
  mt_genes <- read.table(files[i+n],header=FALSE,sep=" ")
  df <- rbind(genes,mt_genes)
  names(df) <- c('chr','start','end','CNAG','score','strand')
  df$Condition <- c(rep('WildType',nrow(genes)),rep('Mutant',nrow(mt_genes)))
  sample <- strsplit(files[i],'_e')[[1]][1]
  title <- paste(strsplit(sample,'_')[[1]][1],strsplit(sample,'_')[[1]][2])
  plotlist[[i]] <- plot_multi_histogram(df, 'score', 'Condition',title)
}

grid.arrange(
  plotlist[[1 ]],plotlist[[8 ]],plotlist[[5 ]],
  plotlist[[3 ]],plotlist[[10]],plotlist[[7 ]],
  plotlist[[2 ]],plotlist[[9 ]],plotlist[[8 ]],
  plotlist[[4 ]],plotlist[[11]],ggplot(),
  ncol=3,nrow=4,
  widths = c(.5,.5,.5),
  top = 'Mutant vs. Wild Type Score Distributions',
  left ="Number of Genes"
  )

#plot_multi_histogram(df, 'score', 'Condition',files[1])
