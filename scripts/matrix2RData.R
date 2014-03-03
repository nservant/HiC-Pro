rm(list=ls())

### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}
require(HiTC)
require(parallel)
options(mc.cores=cpu)

## Import matrices in R
data <- lapply(list.files(matDir, pattern="*.matrix", full.names=TRUE), import.my5C)
data <- HTClist(unlist(data))

## Save raw matrices
assign(obj, data)
save(list=obj, file=file.path(rdataDir, paste(obj, ".RData", sep="")))

## Save annotated matrices
require(rtracklayer)
#source("scripts/hicNorm.R")

if (org == "mm9"){
    genomePack="BSgenome.Mmusculus.UCSC.mm9"
    mapFile="annotation/wgEncodeCrgMapabilityAlign100mer_mm9.bw"
}else if(org == "hg18"){
    genomePack="BSgenome.Hsapiens.UCSC.hg18"
    mapFile="annotation/wgEncodeCrgMapabilityAlign100mer_hg18.bw"
}else if(org == "hg19"){
    genomePack="BSgenome.Hsapiens.UCSC.hg19"
    mapFile="annotation/wgEncodeCrgMapabilityAlign100mer_hg19.bw"
}

map <- import(mapFile, format="BIGWIG", asRangedData=FALSE)
datannot <- setGenomicFeatures(data, genomePack=genomePack, mappability=map)
objannot <- paste(obj, "annot", sep="_")
assign(objannot, datannot)
save(list=objannot, file=file.path(rdataDir, paste(objannot, ".RData", sep="")))
