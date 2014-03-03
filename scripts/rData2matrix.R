rm(list=ls())

### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## cpu
## matDir
## rdata

require(HiTC)
require(parallel)
options(mc.cores=cpu)

## Import matrices in R
load(rdata)
xname <- sub(".RData","", basename(rdata))
x <- eval(as.name(xname))
print(x)

stopifnot(inherits(x,"HTClist"))
mclapply(x, function(xx, org){
    if (isIntraChrom(xx)){
        outname <- paste(xname,".",seqlevels(xx),"_", org,".",seqlevels(xx),"_",org,".matrix",sep="")
    }else{
        outname <- paste(xname,".",seqlevels(xx)[1],"_", org,".",seqlevels(xx)[2],"_",org,".matrix",sep="")
    }
    print(outname)
    export.my5C(x=xx, file=file.path(matDir,outname), genome=org, format='mat', header=TRUE)
}, mc.cores = cpu, org=org)
