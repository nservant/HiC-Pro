## HiC-Pro
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Plot 3C filtering results
##

rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## getHiCMat
## Generate data.frame for ggplot2 graphical output
## x = vector with all expected Hi-C results
##
getHiCMat <- function(x){
  require(RColorBrewer)
  
  invalid.lab <- intersect(names(x), c("Self_Cycle_pairs", "Dangling_end_pairs", "Religation_pairs", "Single-end_pairs", "Dumped_pairs", "Filtered_pairs"))
  valid.lab <- intersect(names(x), c("Valid_interaction_pairs_FF", "Valid_interaction_pairs_RR", "Valid_interaction_pairs_RF", "Valid_interaction_pairs_FR"))
  x <- x[c("Valid_interaction_pairs", valid.lab, invalid.lab)]
  
  ## Add number of invalid pairs
  n.invalid_pairs <-sum(x[invalid.lab])
  x["Invalid_pairs"]<-n.invalid_pairs

  n.valid_pairs <- x["Valid_interaction_pairs"]
  n <- n.invalid_pairs  + n.valid_pairs
  
  ## Get percentage
  x.perc <- round(x/n*100, 1)

  ## multiple plots
  p <- rep(1, length(x))
  p[names(x)%in%valid.lab] <- 2
  p[names(x)%in%invalid.lab] <- 3

  mmat <- data.frame(cbind(lab=names(x), p, count=x, perc=x.perc), stringsAsFactors=FALSE)

  lab.order <- c("Invalid_pairs", "Valid_interaction_pairs", "Valid_interaction_pairs_FF", "Valid_interaction_pairs_RR",
             "Valid_interaction_pairs_RF","Valid_interaction_pairs_FR", "Filtered_pairs", "Dumped_pairs", "Self_Cycle_pairs",
             "Religation_pairs", "Single-end_pairs", "Dangling_end_pairs")
  mmat$lab <- factor(mmat$lab, levels=lab.order)
  #mmat <- mmat[rev(lab.order),]

  ## pos for label
  mmat$pos <- rep(0, length(x))
  for (i in unique(mmat$p)){
      idx <-  which(mmat$p==i)
      counts <- mmat$count[idx]
      names(counts) <- mmat$lab[idx]
      ## get reverse order from levels for position
      s <- counts[rev(intersect(levels(mmat$lab), mmat$lab[idx]))]
      values <- as.numeric(as.character(s))
      mmat$pos[match(names(s), mmat$lab)] <- cumsum(values)-values/2
  }  
  mmat$pos[which(mmat$count==0)] <- NA

  ## Colours
  sel.val <- brewer.pal(8,"Blues")
  sel.inval <- brewer.pal(8,"Reds")

  col <- rep(NA, dim(mmat)[1])
  names(col) <- names(x)
  col[invalid.lab] <- sel.inval[1:length(invalid.lab)]
  col[valid.lab] <- sel.val[1:length(valid.lab)]
  col["Invalid_pairs"] <- "darkgray"
  col["Valid_interaction_pairs"] <- sel.val[5]
  mmat$selcol <- col
      
  #mmat[order(mmat$p), ]
  mmat
}

## plotHiCStat
## Generate ggplot2 plot
## mat = data.frame with all values to plot. see getHiCMat()
## xlab = character for xlabel
## legend = logical. If true, the legend is plotted
##
plotHiCStat <- function(mat, xlab="", legend=TRUE){
  require(RColorBrewer)
  require(ggplot2)
  require(grid)

  ## update labels for plot
  #mat$lab <- paste0(gsub("_", " ", mat$lab)," (%)")

  gp <-ggplot(mat, aes(x=p, as.numeric(count), fill=lab)) + geom_bar(width=.7,stat="identity", colour="gray") + theme_minimal() + 
      theme(axis.title=element_text(face="bold", size=6), axis.ticks = element_blank(),  axis.text.y = element_text(size=5), axis.text.x = element_text(size=6))+
        xlab(xlab) + ylab("Read Counts") +
            scale_x_discrete(breaks=c("1", "2", "3"), labels=c("All Pairs","Valid 3C Pairs","Invalid 3C Pairs"))+
              geom_text(aes(x=p, y=as.numeric(pos), label=paste0(perc,"%")),fontface="bold", size=2) +
                ggtitle("Statistics of Read Pairs Alignment on Restriction Fragments") + theme(plot.title = element_text(lineheight=.8, face="bold", size=6))

  if (legend){
    scol <- mat$selcol
    names(scol) <- mat$lab
    gp = gp + scale_fill_manual(values=scol) + guides(fill=guide_legend(title="")) + theme(plot.margin=unit(x=c(1,0,0,0), units="cm"), legend.position="right",
                                                        legend.margin=margin(.5,unit="cm"), legend.text=element_text(size=5))
  }else{
    gp = gp + scale_fill_manual(values=as.character(col)) + theme(plot.margin=unit(c(1,0,1.9,0),"cm"))+ guides(fill=FALSE)
  }
  gp
}

####################################
##
## plotMappingPortion.R
##
####################################

## Get HiC stat files for all fastq files of a given sample
allrsstat <- list.files(path=hicDir, pattern=paste0("^[[:print:]]*\\.RSstat$"), full.names=TRUE)
print(hicDir)
print(allrsstat)
stopifnot(length(allrsstat)>0)

## Get statistics summary
stats_per_fastq<- lapply(allrsstat, read.csv, sep="\t", as.is=TRUE, comment.char="#", header=FALSE, row.names=1)
stats_per_sample<- rowSums(do.call(cbind,stats_per_fastq))
print (stats_per_sample)

## Make plots
mat <- getHiCMat(stats_per_sample)
p1 <- plotHiCStat(mat, xlab=sampleName)
ggsave(filename=file.path(picDir, paste0("plotHiCFragment_",sampleName,".pdf")), p1, width=5, height=5)
