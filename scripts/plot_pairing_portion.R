## HiC-Pro
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Plot R1/R2 Pairing results
##

rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## getPairMat
## Generate data.frame for ggplot2 graphical output
## x = vector with all expected Hi-C results
##
getPairMat <- function(x, x.perc, rmMulti=0, rmSingle=0){
  require(RColorBrewer)

  n.plots <- 4
  notreported.lab <- c("Low_qual_pairs")
  ## Remove multi pairs
  if (rmMulti == 1 && rmSingle ==1){
    n.plots <- n.plots + 2
    notreported.lab <- c(notreported.lab, "Multiple_pairs_alignments", "Pairs_with_singleton")
  }else if (rmMulti == 1 && rmSingle == 0){
      n.plots <- n.plots + 2
      notreported.lab <- c(notreported.lab, "Multiple_pairs_alignments", "Multiple_singleton_alignments")
  }else if (rmMulti == 0 && rmSingle == 1){
      n.plots <- n.plots + 1
      notreported.lab <- c(notreported.lab, "Pairs_with_singleton")
  }
  
  reported.lab <- "Reported_pairs"
  allnotreported.lab <- "Not_Reported_pairs"
  un.lab <- "Unmapped_pairs"
  
  n.reported_pairs <- x[reported.lab]
  n.notreported_pairs <-sum(x[notreported.lab])

  ## Add total number of not reported pairs
  x[allnotreported.lab]<-n.notreported_pairs
  x.perc[allnotreported.lab] <- round(n.notreported_pairs/x["Total_pairs_processed"]*100,1)

  ## Check
  print (x[un.lab])
  print (x[reported.lab])
  print(x[allnotreported.lab])
  stopifnot(x[un.lab]+x[reported.lab]+x[allnotreported.lab]==x["Total_pairs_processed"])
    
  ## Get percentage
  x.perc <- round(x.perc, 1)
  
  ## multiple plots
  p <- rep(1, n.plots)
  print(p)
  print(c(reported.lab, allnotreported.lab, un.lab, notreported.lab))
  names(p) <- c(reported.lab, allnotreported.lab, un.lab, notreported.lab)
  p[notreported.lab] <- 2

  mmat <- data.frame(cbind(lab=names(p), p, count=x[names(p)], perc=x.perc[names(p)]), stringsAsFactors=FALSE)

  ## pos for label
  mmat$pos <- rep(0, length(p))
  for (i in unique(mmat$p)){
    idx <-  which(mmat$p==i)
    mmat$pos[idx] <- cumsum(as.numeric(as.character(mmat$count[idx])))-as.numeric(as.character(mmat$count[idx]))/2
  }  
  mmat$pos[which(mmat$count==0)] <- NA

  ## Colours
  sel.reported <- brewer.pal(6,"Blues")
  sel.notreported <- brewer.pal(6,"Reds")

  col <- rep(NA, dim(mmat)[1])
  names(col) <- names(p)
  col[notreported.lab] <- sel.notreported[1:length(notreported.lab)]
  col[c(reported.lab, allnotreported.lab)] <- sel.reported[1:2]
  col[un.lab] <- "darkgray"
  mmat$selcol <- col
  
  mmat[order(mmat$p), ]
}

## plotPairStat
## Generate ggplot2 plot
## mat = data.frame with all values to plot. see getHiCMat()
## xlab = character for xlabel
## legend = logical. If true, the legend is plotted
##
plotPairStat <- function(mat, xlab="", legend=TRUE){
    require(RColorBrewer)
    require(ggplot2)
    require(grid)
    
    gp <-ggplot(mat, aes(x=p, as.numeric(count), fill=lab)) + geom_bar(width=.7,stat="identity", colour="gray") + theme_minimal() +
            theme(axis.title=element_text(face="bold", size=6), axis.ticks = element_blank(),  axis.text.y = element_text(size=5), axis.text.x = element_text(size=6))+
                xlab(xlab) + ylab("Read Counts") +
                    scale_x_discrete(breaks=c("1", "2"), labels=c("All Pairs","Filtered Pairs"))+
                        geom_text(aes(x=p, y=as.numeric(pos), label=paste(perc,"%")),fontface="bold", size=2) +
                            ggtitle("Statistics after read pairing") + theme(plot.title = element_text(lineheight=.8, face="bold", size=6))
    
    if (legend){
        scol <- mat$selcol
        names(scol) <- mat$lab
        gp = gp + scale_fill_manual(values=scol) + guides(fill=guide_legend(title="")) + theme(plot.margin=unit(x=c(1,0,0,0), units="cm"),
                                                              legend.position="right", legend.margin=margin(.5, unit="cm"),legend.text=element_text(size=5))
    }else{
        gp = gp + scale_fill_manual(values=as.character(col)) + theme(plot.margin=unit(c(1,0,1.9,0),"cm"))+ guides(fill=FALSE)
    }
    gp
}

####################################
##
## plotPairingPortion.R
##
####################################

## Get HiC stat files for all fastq files of a given sample
allpairstat <- list.files(path=bwtDir, pattern=paste0("^[[:print:]]*\\.pairstat$"), full.names=TRUE)
stopifnot(length(allpairstat)>0)

## Get statistics summary
stats_per_fastq<- lapply(allpairstat, read.csv, sep="\t", as.is=TRUE, header=FALSE, comment.char="#", row.names=1)
stats_per_sample<- rowSums(do.call(cbind,lapply(stats_per_fastq, "[", 1)))
perc_per_sample<- rowMeans(do.call(cbind,lapply(stats_per_fastq, "[", 2)))

print(stats_per_sample)
## Make plots
mat <- getPairMat(stats_per_sample, perc_per_sample, rmMulti=rmMulti, rmSingle=rmSingle)
mat$lab <- factor(mat$lab, levels=c("Unmapped_pairs", "Not_Reported_pairs", "Reported_pairs", "Low_qual_pairs",  "Pairs_with_singleton", "Multiple_pairs_alignments"))
p1 <- plotPairStat(mat, xlab=sampleName)
ggsave(file.path(picDir, paste0("plotMappingPairing_",sampleName,".pdf")), p1, width=5, height=5)
