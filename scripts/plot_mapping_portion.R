## Nicolas Servant
## HiC-Pro
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Plot mapping proportion
##
rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## Multiple plot function
##
## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
##
## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.
##
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
      print(plots[[1]])
      
  } else {
      ## Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## getMapMat
## Generate data.frame for ggplot2 graphical output
## n = numeric. Total number of reads
## n.map = numeric. Total number of aligned reads
## n.map.glob = numeric. Total number of aligned reads with global mapping
## n.map.loc = numeric. Total number of aligned reads with local mapping
##
getMapMat <- function(n, n.map, n.map.glob, n.map.loc){

  ## check
  stopifnot(n.map.glob+n.map.loc==n.map)
  n.unmap = n - n.map
  
  ## Get percentage
  n.map.perc <- round(100*n.map/n, digits=1)
  n.unmap.perc <- round(100*n.unmap/n, digits=1)
  n.map.glob.perc <-  round(100*n.map.glob/n, digits=1)
  n.map.loc.perc <-  round(100*n.map.loc/n, digits=1)

  p <- c(rep("1",2), rep("2",3))
  count <- c(n.map, n.unmap, n.map.glob, n.map.loc, n.unmap)
  perc <- c(n.map.perc, n.unmap.perc, n.map.glob.perc, n.map.loc.perc, n.unmap.perc)
  lab <- c("n.map","n.un","n.glob","n.loc","n.un")
  mmat <- data.frame(cbind(lab, p, count, perc), stringsAsFactors=FALSE)
  mmat$pos <- as.vector(unlist(sapply(unique(mmat$p), function(i){
    idx <-  which(mmat$p==i)
    cumsum(as.numeric(as.character(mmat$count[idx])))-as.numeric(as.character(mmat$count[idx]))/2
  })))

  mmat
}

## ploMapStat
## Generate ggplot2 plot
## mat = data.frame for ggplot2 input. see getMapMat()
## xlab = character for xlabel
## legend = logical. If true, the legend is plotted
##
ploMapStat <- function(mat, sampleName="", tag="", legend=TRUE){
  require(RColorBrewer)
  require(ggplot2)
  require(grid)
  
  sel.colours <- brewer.pal(6,"Blues")
  tit <- "Statistics of Read Alignments"
  if (tag != ""){
    tit <- paste0(tit," - ", tag," Tags") 
  }
   
  
  gp <- ggplot(mat, aes(x=p, as.numeric(count), fill=lab)) + geom_bar(width=.7,stat="identity", colour="gray") + theme_minimal() + 
          theme(axis.title=element_text(face="bold", size=6), axis.ticks = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank()) +
              xlab(sampleName) + ylab("Read Counts")  +
                  geom_text(aes(x=p, y=as.numeric(pos), label=paste(perc,"%")),fontface="bold", size=2)+
                      ggtitle(tit) + theme(plot.title = element_text(lineheight=.8, face="bold", size=8))

  if (legend){
    gp = gp + scale_fill_manual(values=c("darkgray", sel.colours[2:4]), labels = c("Not aligned (%)", "Trimmed read Mapping (%)","Full read mapping (%)","Aligned reads (%)")) + guides(fill=guide_legend(title="")) + theme(plot.margin=unit(x=c(1,0,0,0), units="cm"), legend.position="bottom", legend.text=element_text(size=5))
  }else{
    gp = gp + scale_fill_manual(values=c("darkgray", sel.colours[2:4])) + theme(plot.margin=unit(c(1,0,1.45,0),"cm"))+ guides(fill=FALSE)
  }
  gp
}                   


####################################
##
## plotMappingPortion.R
##
####################################
## Get Mapping stat files for all fastq files of a given sample
allmapstat_r1 <- list.files(path=bwtDir, pattern=paste0("^[[:print:]]*",r1tag,"[[:print:]]*\\.mapstat$"), full.names=TRUE)
allmapstat_r2 <- list.files(path=bwtDir, pattern=paste0("^[[:print:]]*",r2tag,"[[:print:]]*\\.mapstat$"), full.names=TRUE)
stopifnot(length(allmapstat_r1)>0 && length(allmapstat_r2)>0)
print(allmapstat_r1)
print(allmapstat_r2)

## Get statistics summary
stats_per_fastq_r1<- sapply(allmapstat_r1, read.csv, sep="\t", row.names=1, as.is=TRUE, comment.char="#", header=FALSE)
stats_per_sample_r1<- colSums(do.call(rbind, stats_per_fastq_r1))
stats_per_fastq_r2<- sapply(allmapstat_r2, read.csv, sep="\t", row.names=1, as.is=TRUE, comment.char="#", header=FALSE)
stats_per_sample_r2<- colSums(do.call(rbind, stats_per_fastq_r2))

print(stats_per_sample_r1)
print(stats_per_sample_r2)

## Make plots
mat_r1 <- getMapMat(stats_per_sample_r1[1], stats_per_sample_r1[2], stats_per_sample_r1[3], stats_per_sample_r1[4])
mat_r1$lab <- factor(mat_r1$lab, levels=c("n.un","n.loc","n.glob","n.map"))
p1 <- ploMapStat(mat_r1, sampleName=sampleName, tag="R1")
mat_r2 <- getMapMat(stats_per_sample_r2[1], stats_per_sample_r2[2], stats_per_sample_r2[3], stats_per_sample_r2[4])
mat_r2$lab <- factor(mat_r2$lab, levels=c("n.un","n.loc","n.glob","n.map"))
p2 <- ploMapStat(mat_r2, sampleName=sampleName, tag="R2", legend=FALSE)
pdf(file.path(picDir, paste0("plotMapping_",sampleName,".pdf")), width=7, height=4)
multiplot(p1, p2, cols=2)
dev.off()

