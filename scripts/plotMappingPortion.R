## Nicolas Servant
## Plot mapping proportion


rm(list=ls())
require(RColorBrewer)

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

bwt_glob <- file.path(bwtDir, "bwt2_global", sampleName)
bwt_loc <- file.path(bwtDir, "bwt2_local", sampleName)
bwt_merge <- file.path(bwtDir, "bwt2", sampleName)

plotMappingStat <- function(bg, bl, bm, xlab, legend=FALSE){
    n <- sum(bm[,2])
    n.map <- sum(bm[which(bm[,1]=="Uniquely_Mapped"),2])
    n.unmap <- n - n.map
    n.map.perc <- round(100*n.map/n, digits=1)
    ng.map <- sum(bg[which(bg[,1]=="Uniquely_Mapped"),2])
    ng.map.perc <-  round(100*ng.map/(n), digits=1)
    nl.map <- sum(bl[which(bl[,1]=="Uniquely_Mapped"),2])
    nl.map.perc <-  round(100*nl.map/(n), digits=1)
    
    stopifnot(ng.map+nl.map == n.map)

    sel.colours <- brewer.pal(6,"Blues")
    mat <- matrix(c(n.map, 0, 0, n.unmap, 0, ng.map, nl.map, n.unmap), ncol=2)
    xpos <- barplot(mat, col=c(sel.colours[6:4],"gray"), ylab="Number of Reads",
                    las=2, xlab=xlab, cex.axis=.7, cex.lab=.7, cex.names=.8, xlim=c(0,4))
    
    mtext(side=1, line=-1.3, at=xpos, text=paste(c(n.map.perc, ng.map.perc),"%",sep=""), font=2)
    text(x=xpos[2],y=ng.map, labels=paste(nl.map.perc,"%",sep=""), font=2, pos=3)
    text(x=xpos,y=n.map, labels=paste(100-n.map.perc,"%",sep=""), font=2, pos=3)

    legend(x="topright", legend=c("Unique matches(%)", "Global mapping(%)", "Local Mapping(%)", "Not matches(%)"), 
           fill=c(sel.colours[6:4], "gray"), cex=0.7, inset=0.01)   
}

writeSummaryTable <- function(bg, bl, bm, tags, ...){
    n <- sum(bm[,2])
    n.map <- sum(bm[which(bm[,1]=="Uniquely_Mapped"),2])
    n.unmap <- n - n.map
    n.map.perc <- round(100*n.map/n, digits=1)
    ng.map <- sum(bg[which(bg[,1]=="Uniquely_Mapped"),2])
    ng.map.perc <-  round(100*ng.map/(n), digits=1)
    nl.map <- sum(bl[which(bl[,1]=="Uniquely_Mapped"),2])
    nl.map.perc <-  round(100*nl.map/(n), digits=1)
    
    out <- as.data.frame(matrix(c(paste("nbreads_", tags, sep=""), n,
                                  paste("nbreads_mapped_", tags, sep=""), n.map,
                                  paste("nbreads_unmapped_", tags, sep=""), n.unmap,
                                  paste("nbreads_mapped_global_", tags, sep=""),ng.map,
                                  paste("nbreads_mapped_local_", tags, sep=""),nl.map)
                                , ncol=2, byrow=TRUE))

    write.table(out, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, ...)
}

## R1 tag
bg.1 <- read.table(file.path(bwt_glob, paste(sampleName,"_R1.mapstat", sep="")))
bl.1 <- read.table(file.path(bwt_loc, paste(sampleName,"_R1.mapstat", sep="")))
bm.1 <- read.table(file.path(bwt_merge, paste(sampleName,"_R1.mapstat", sep="")))

## R2 tag
bg.2 <- read.table(file.path(bwt_glob, paste(sampleName,"_R2.mapstat", sep="")))
bl.2 <- read.table(file.path(bwt_loc, paste(sampleName,"_R2.mapstat", sep="")))
bm.2 <- read.table(file.path(bwt_merge, paste(sampleName,"_R2.mapstat", sep="")))

png(file.path(picDir,"plotGenomeMapping.png"), units="in", res=300, heigh=5, width=9)
par(font.lab=2, mai=c(1.2, 1.1, 0.2, 0.1), mfrow=c(1,2))
plotMappingStat(bg.1, bl.1, bm.1, xlab=paste(sampleName,"- R1 Tags"))
plotMappingStat(bg.2, bl.2, bm.2, xlab=paste(sampleName,"- R2 Tags"))
dev.off()

writeSummaryTable(bg.1, bl.1, bm.1, tags="R1", file=file.path(docDir,paste(sampleName,"mappingstat.tsv", sep="_")))
writeSummaryTable(bg.2, bl.2, bm.2, tags="R2", file=file.path(docDir,paste(sampleName,"mappingstat.tsv", sep="_")), append=TRUE)
