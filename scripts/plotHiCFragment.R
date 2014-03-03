rm(list=ls())
require(RColorBrewer)

### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

plotFragmentInformation <- function(headFile){

    fragmentInfo <- read.table(headFile, row.names=1)
    fragmentInfo[,1]<-as.numeric(gsub(",","",fragmentInfo[,1]))
    
    n <- fragmentInfo["bothSideMapped",1]
    n.val <- fragmentInfo["validPairs",1]
    n.val.perc <- round(100*n.val/n, digits=1)
    n.inval <- fragmentInfo["invalidPairs",1]
    n.inval.perc <- round(100*n.inval/n, digits=1)
    n.error <- fragmentInfo["errorPairs",1]
    n.error.perc <- round(100*n.error/n, digits=1)
    
    ## Details of valid pairs
    n.val.1 <- fragmentInfo["different|->.->",1]
    n.val.1.perc <- round(100*n.val.1/n.val, digits=1)
    n.val.2 <- fragmentInfo["different|->.<-",1]
    n.val.2.perc <- round(100*n.val.2/n.val, digits=1)
    n.val.3 <- fragmentInfo["different|<-.->",1]
    n.val.3.perc <- round(100*n.val.3/n.val, digits=1)
    n.val.4 <- fragmentInfo["different|<-.<-",1]
    n.val.4.perc <- round(100*n.val.4/n.val, digits=1)
    
    ## Details of invalid pairs
    n.inval.2 <- fragmentInfo["selfCircle",1]
    n.inval.2.perc <- round(100*n.inval.2/n.inval, digits=1)
    n.inval.1 <- fragmentInfo["danglingEnd",1]
    n.inval.1.perc <- round(100*n.inval.1/n.inval, digits=1)
    
    sel.colours.b <- brewer.pal(9,"Blues")
    sel.colours.r <- brewer.pal(9,"Reds")
    sel.colours <- c(sel.colours.b[8], sel.colours.r[8],"darkgreen", sel.colours.b[6:3], sel.colours.r[6:5])
    
    mat <- matrix(c(n.val, n.inval, n.error,0,0,0,0,0,0,
                    0, 0, 0, n.val.1, n.val.2, n.val.3, n.val.4,0,0,
                    0, 0, 0, 0, 0, 0, 0, n.inval.1, n.inval.2), ncol=3)

    xpos <- barplot(mat, col=sel.colours, ylab="Number of Reads",
                    las=1, xlab=sampleName, names.arg=c("allPairs","validPairs", "invalidPairs"), cex.axis=.7, cex.lab=.7, cex.names=.8, xlim=c(0,4))
    
    mtext(side=1, line=-1.3, at=xpos, text=paste(c(n.val.perc, n.val.1.perc, n.inval.1.perc),"%",sep=""), font=2)
    
    text(x=xpos[1],y=n.val, labels=paste(n.inval.perc,"%",sep=""), font=2, pos=3)
    
    text(x=xpos[2],y=n.val.1, labels=paste(n.val.2.perc,"%",sep=""), font=2, pos=3)
    text(x=xpos[2],y=n.val.1+n.val.2, labels=paste(n.val.3.perc,"%",sep=""), font=2, pos=3)
    text(x=xpos[2],y=n.val.1+n.val.2+n.val.3, labels=paste(n.val.4.perc,"%",sep=""), font=2, pos=3)
    
    text(x=xpos[3],y=n.inval.1, labels=paste(n.inval.2.perc,"%",sep=""), font=2, pos=3)
    
    legend(x="topright",
           legend=c("validPairs(%)", "invalidPairs(%)", "errorPairs", "validPairs(->->)(%)", "validPairs(-><-)(%)", "validPairs(<-->)(%)", "validPairs(<-<-)(%)", "danglingEnd(%)", "selfCircle(%)"), 
           fill=sel.colours, cex=0.7, inset=0.01)   
}

headFile <- list.files(dataDir, pattern="*.head", full.names=TRUE)
png(file.path(picDir,"plotFragmentInfo.png"), units="in", res=300, heigh=5, width=6)
par(font.lab=2, mai=c(1.2, 1.1, 0.2, 0.1))
plotFragmentInformation(headFile[1])
dev.off()

