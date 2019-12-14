library(dplyr)

enrich_files <- list.files("~/Projects/fixpeaks/sums3-enriched/",full.names = TRUE)
all_g_files  <- list.files("~/Projects/fixpeaks/sums-peak-by-gene/",full.names = TRUE)

# plotlist     <- list()
# percent_fold <- list()
plot_scatters <- function(df,title,subtitle1,subtitle2){   
  sc <- ggplot(df, aes(x=wt_score,y=mut_score)) + 
    geom_point(alpha = 0.18, aes(x=wt_score, y=mut_score, color=color), show.legend = FALSE) +
    geom_abline(slope=1,intercept=0,col="grey") +
    scale_x_log10( limits=c(339,65905)) + 
    scale_y_log10(limits=c(339,65905)) +
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

scatterlist <- list()
all  <- data.frame(a=NA,b=NA,c=NA,CNAG=NA,wt_score=NA,g=NA,d=NA,e=NA,f=NA,mut_score=NA,h=NA,combined=NA,color=NA,name=NA)
diff <- data.frame(CNAG=NA, wt_score=NA, mut_score=NA, combined=NA, color=NA, name=NA)

n    <- length(enrich_files)/2


for (i in 1:n){
  enriched_wt         <- read.table(enrich_files[i],header=FALSE,sep=" ")
  enriched_mut        <- read.table(enrich_files[i+n],header=FALSE,sep=" ")
  all_wt         <- read.table(all_g_files[i],header=FALSE,sep=" ")
  all_mut        <- read.table(all_g_files[i+n],header=FALSE,sep=" ")
  
  enriched_wt[,c(1,2,3,6)]  <- NULL
  enriched_mut[,c(1,2,3,6)] <- NULL
  all_wt[,c(1,2,3,6)]  <- NULL
  all_mut[,c(1,2,3,6)] <- NULL
  
  names(enriched_wt)  <- c('CNAG','wt_score')
  names(enriched_mut) <- c('CNAG','mut_score')
  names(all_wt)  <- c('CNAG','wt_score')
  names(all_mut) <- c('CNAG','mut_score')
  nrow(enriched_wt)
  nrow(enriched_mut)
  nrow(all_wt)
  nrow(all_mut)
  
  # get the mutual gene list and grab scores 
  enriched_both <- full_join(enriched_wt,enriched_mut,by='CNAG')
  enriched_both <- left_join(enriched_both,all_wt,by='CNAG')
  enriched_both <- left_join(enriched_both,all_mut,by='CNAG')
  
  enriched_both[,c(2,3)]       <- NULL
  names(enriched_both)[c(2,3)] <- c('wt_score','mut_score')

  enriched_both$combined <- enriched_both$wt_score/enriched_both$mut_score
  enriched_both$color    <- 'darkblue'
  enriched_both$color[ enriched_both$combined>1.5 | enriched_both$combined<.67 ] <- 'red'
  
  title        <- strsplit(strsplit(enrich_files[i],"//")[[1]][2], '_e')[[1]][1]
  total        <- length(enriched_both$combined)
  
  num          <- sum(enriched_both$combined>1.5)
  percent_fold <- round(num/total,4)
  subtitle     <- paste(paste(num,"/",total,sep=""),"genes: FC > 1.5 (",percent_fold,")")
  
  num2          <- sum(enriched_both$combined<.67)
  percent_fold2 <- round(num2/total,4)
  subtitle2     <- paste(paste(num2,"/",total,sep=""),"genes: FC < .67 (",percent_fold2,")")
  
  scatterlist[[i]] <- plot_scatters(enriched_both,title,subtitle,subtitle2)
  
  different <- enriched_both[enriched_both$combined>1.5,]
  different$name <- title
  diff <- rbind(diff,different)
  joined$name <- title
  all <- rbind(all,joined)
  print(max(enriched_both$mut_score))
  print(max(enriched_both$wt_score)) 
}

grid.arrange(
  scatterlist[[1 ]],scatterlist[[9 ]],scatterlist[[5 ]],
  scatterlist[[3 ]],scatterlist[[11]],scatterlist[[7 ]],
  scatterlist[[2 ]],scatterlist[[10 ]],scatterlist[[6 ]],
  scatterlist[[4 ]],scatterlist[[12]],scatterlist[[8 ]],
  ncol=3,nrow=4,
  widths = c(1,1,1),
  top = 'WT vs. Mutant Scores (Using Sums & k=3, Joining Enriched)',
  left ="Mutanat Score",
  bottom="Wild Type Score"
)
