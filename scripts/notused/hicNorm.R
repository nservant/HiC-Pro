### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
      eval(parse(text=args[[i]]))
}

## CPU
require(parallel)
require(HiTC)

## Load RData file
print(rdata)
load(rdata)
xname <- sub(".RData","", basename(rdata))
x <- eval(as.name(xname))
print(x)

stopifnot(inherits(x,"HTClist"))
nbchrom <- length(x)

## Normalisation
if (norm == "ICE"){
    objout <- paste(xname,"_iced",sep="")
    z <- mclapply(x, normICE, mc.cores = cpu)
    stopifnot(length(z)==nbchrom)
    z <- HTClist(z)
    assign(objout, z)
    save(list=objout, file=file.path(outDir,paste(xname,"_iced.RData", sep="")))
}else if (norm == "LGF"){
    stopifnot(length(intersect(colnames(elementMetadata(x_intervals(x$chr1chr1))), c("len","GC")))==2)
    objout <- paste(xname,"_lgf",sep="")
    z <- mclapply(x, normLGF, mc.cores = cpu, family="nb")
    stopifnot(length(z)==nbchrom)
    z <- HTClist(z)
    assign(objout, z)
    save(list=objout, file=file.path(outDir,paste(xname,"_lgf.RData", sep="")))
}
print(z)

