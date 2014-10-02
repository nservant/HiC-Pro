rm(list=ls())

### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

matrix.import <- function(filename=""){
  if(!file.exists(filename)){
    error.msg <- paste("can't fild matrix file -- ",filename ,sep="")
    stop(error.msg)
  }
  cur.matrix <- try(read.table(filename))
  if(class(cur.matrix)=='try-error'){
    return (0)
  }
  else{
    return (as.matrix(cur.matrix))
  }
}


perfect.cis.filenameA <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genome1,".matrix",sep="")
perfect.cis.filenameB <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genome2,".matrix",sep="")
perfect.trans.filename <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genome2,".matrix",sep="")

common.cis.filenameA <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genomeCommon,".cis.matrix",sep="")
common.cis.filenameB <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genomeCommon,".cis.matrix",sep="")
common.trans.filenameA <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genomeCommon,".trans.matrix",sep="")
common.trans.filenameB <- paste(commonDir,"/",sampleName,".",chrmstr2,"_",genome2,".",chrmstr1,"_",genomeCommon,".trans.matrix",sep="")

perfect.cisA <- matrix.import(filename=perfect.cis.filenameA)
perfect.cisB <- matrix.import(filename=perfect.cis.filenameB)
perfect.trans <- matrix.import(filename=perfect.trans.filename)
common.cisA <- matrix.import(filename=common.cis.filenameA)
common.cisB <- matrix.import(filename=common.cis.filenameB)
common.transA <- matrix.import(filename=common.trans.filenameA)
common.transB <- matrix.import(filename=common.trans.filenameB)

common.transB <- t(common.transB)

all.cisA <- perfect.cisA + common.cisA
all.cisB <- perfect.cisB + common.cisB
all.trans <- perfect.trans + common.transA + common.transB

cisA.output.name <- paste(outDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genome1,".cis.matrix",sep="")
cisB.output.name <- paste(outDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genome2,".cis.matrix",sep="")
trans.output.name <- paste(outDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genome2,".trans.matrix",sep="")

write.table(all.cisA,file=cisA.output.name,quote=FALSE,sep="\t")
write.table(all.cisB,file=cisB.output.name,quote=FALSE,sep="\t")
write.table(all.trans,file=trans.output.name,quote=FALSE,sep="\t")
