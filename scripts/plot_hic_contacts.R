## HiC-Pro
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Plot Valid Interaction Statistics
##

rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## getContactsStatMat
## Generate data.frame for ggplot2 graphical output
## x = vector with all expected Hi-C results
##

getContactsStatMat <- function(x){
  require(RColorBrewer)

  ## counts  
  ncontacts <- x["valid_interaction",1]
  rmdup <- x["valid_interaction_rmdup",1]
  ndup <- x["valid_interaction",1] - x["valid_interaction_rmdup",1]
  stopifnot(rmdup+ndup==ncontacts)

  cis_inter <- x["cis_interaction",1]
  cis_sr <- x["cis_shortRange",1]
  cis_lr <-  x["cis_longRange",1]
  trans_inter <- x["trans_interaction",1]
  stopifnot(sum(c(cis_sr,cis_lr), na.rm=TRUE)==cis_inter)

  count <- c(rmdup, ndup, trans_inter, cis_sr, cis_lr)

  ## perc
  rmdup.perc <- round(100*rmdup/ncontacts)
  ndup.perc <- round(100*ndup/ncontacts)
  trans_inter.perc <- round(100*trans_inter/ncontacts)
  cis_sr.perc <- round(100*cis_sr/ncontacts)
  cis_lr.perc <- round(100*cis_lr/ncontacts)
  
  perc <- c(rmdup.perc, ndup.perc, trans_inter.perc, cis_sr.perc, cis_lr.perc)
  p <- c(rep("1",2), rep("2", 3))

  lab <- c("n.a.rmdup", "n.b.dup", "n.c.trans", "n.d.cis.sr", "n.e.cis.lr")
  mmat <- data.frame(cbind(lab, p, count, perc), stringsAsFactors=FALSE)
  mmat$pos <- as.vector(unlist(sapply(unique(mmat$p), function(i){
      idx <-  which(mmat$p==i)
      values <- as.numeric(as.character(mmat$count[idx]))
      cumsum(values)-values/2
  })))  

  mmat$lab <- factor(mmat$lab, levels=c("n.b.dup", "n.a.rmdup", "n.e.cis.lr", "n.d.cis.sr", "n.c.trans"))
  mmat
}

## plotDedup
## Generate ggplot2 plot
## mat = data.frame for ggplot2 input. see getContactsStatMat()
## xlab = character for xlabel
## legend = logical. If true, the legend is plotted
##
plotDedup <- function(mat, sampleName="", legend=TRUE){
  require(RColorBrewer)
  require(ggplot2)
  require(grid)
  
  sel.colours <- brewer.pal(12,"Paired")[c(8,7,2,1,3)] 

  gp <- ggplot(mat, aes(x=p, as.numeric(count), fill=lab)) +
    geom_bar(width=.7,stat="identity", colour="gray") + theme_minimal() + 
      theme(axis.title=element_text(face="bold", size=6), axis.ticks = element_blank(), axis.text.y = element_text(size=5), axis.text.x = element_blank()) +
          xlab(sampleName) + ylab("Read Counts")  +
            geom_text(aes(x=p, y=as.numeric(pos), label=paste(perc,"%")),fontface="bold", size=2) +
                ggtitle("Valid Pairs - duplicates and contact ranges") + theme(plot.title = element_text(lineheight=.8, face="bold", size=6))

  if (legend){
    gp = gp + scale_fill_manual(values=sel.colours, labels = c("Duplicates (%)", "Valid Interactions (%)", "Cis long-range (>20kb) (%)", "Cis short-range contacts (<20kb) (%)", "Trans Contacts (%)")) +
                                    guides(fill=guide_legend(title="")) +
                                        theme(plot.margin=unit(x=c(1,0,0,0), units="cm"), legend.position="right", legend.margin=margin(.5,unit="cm"),
                                              legend.text=element_text(size=5))
  }else{
    gp = gp + scale_fill_manual(values=sel.colours) + theme(plot.margin=unit(c(1,0,1.9,0),"cm"))+ guides(fill=FALSE)
  }
  gp
}

plotDistanceHist <- function(mat, sampleName="", n=""){
  require(RColorBrewer)
  require(ggplot2)
  require(grid)
  
  gp <- ggplot(mat, aes(x=mids, y=allcounts))+
     geom_bar(stat="identity", alpha=.5, color="darkgray", fill="blue4")+theme_minimal()+
      theme(axis.title=element_text(face="bold", size=6), axis.text.y = element_text(size=5), axis.text.x = element_text(size=5)) + 
      scale_x_continuous(breaks=c(seq(0, 500, by=50), seq(from = 600, to = 1500, by = 200), 1500), labels=c(seq(0, 500, by=50), seq(from = 600, to = 1500, by = 200), ">1500"))+
      xlab(sampleName) + ylab(paste0("Read Counts - subset of ", n, " interactions")) +
      ggtitle("Valid Pairs - Fragment size distribution") + theme(plot.title = element_text(lineheight=.8, face="bold", size=6))
gp
}


####################################
##
## plot_hic_contacts.R
##
####################################

## Get HiC stat files for all fastq files of a given sample
mergestat <- list.files(path=statsDir, pattern=paste0("^[[:print:]]*\\.mergestat$"), full.names=TRUE)
print(mergestat)
stopifnot(length(mergestat)>0)

## Get statistics summary
stats_per_sample<- read.csv(mergestat, sep="\t", as.is=TRUE, comment.char="#", header=FALSE, row.names=1)
print(stats_per_sample)
mat <- getContactsStatMat(stats_per_sample)
p1 <- plotDedup(mat, sampleName)
ggsave(filename=file.path(picDir, paste0("plotHiCContactRanges_",sampleName,".pdf")), p1, width=5, height=5)


## Histogram of insert size
allvalidpairs <- list.files(path=hicDir, pattern=paste0("^[[:print:]]*\\.validPairs$"), full.names=TRUE)
stats_per_validpairs<- lapply(allvalidpairs, read.csv, sep="\t", as.is=TRUE, header=FALSE, row.names=1, nrow=100000)
lv <- sapply(stats_per_validpairs, "[", 7)
lv <- lapply(lv, function(x){as.numeric(x[which(x!="None" & ! is.na(x))])})
allhist <- lapply(lv, hist, breaks=c(seq.int(from=0, to=1500, by=10), Inf), plot=FALSE)
allcounts <- Reduce("+", lapply(allhist, "[[", "counts"))

if (max(allcounts)>0){
  mids <- allhist[[1]]$mids
  mat<-data.frame(allcounts=allcounts, mids=mids)
  mat[dim(mat)[1],2]<-1505
  print(allcounts)
  p2 <- plotDistanceHist(mat, sampleName, n=100000*length(allvalidpairs))
  ggsave(filename=file.path(picDir, paste0("plotHiCFragmentSize_",sampleName,".pdf")), p2, width=7, height=5)
}

