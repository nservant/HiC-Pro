rm(list=ls())

### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

matrix.import <- function(filename="",type="prob"){
  if(!file.exists(filename)){
    error.msg <- paste("can't fild matrix file -- ",filename ,sep="")
    stop(error.msg)
  }
  cur.matrix <- try(read.table(filename))
  if(class(cur.matrix)=='try-error'){
    if(type=="prob"){
      return (1)
    }
    else{
      return (0)
    }
  }
  else{
    return (as.matrix(cur.matrix))
  }
}

cis.prob.filename <- paste(probDir,"/",sampleName,".",chrmstr1,".",chrmstr2,".cis.prob.matrix",sep="")
common.filenameA <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genomeCommon,".matrix",sep="")
common.filenameB <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genomeCommon,".matrix",sep="")

cis.prob.matrix <- matrix.import(filename=cis.prob.filename,type="prob")
trans.prob.matrix <- 1-cis.prob.matrix

common.matrixA <- matrix.import(filename=common.filenameA,type="other")
common.matrixB <- matrix.import(filename=common.filenameB,type="other")

common.matrixA.cis <- common.matrixA*cis.prob.matrix
common.matrixB.cis <- common.matrixB*cis.prob.matrix

common.matrixA.trans <- common.matrixA*trans.prob.matrix
common.matrixB.trans <- common.matrixB*trans.prob.matrix

cisA.output.name <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genomeCommon,".cis.matrix",sep="")
cisB.output.name <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genomeCommon,".cis.matrix",sep="")
transA.output.name <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genomeCommon,".trans.matrix",sep="")
transB.output.name <- paste(commonDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genomeCommon,".trans.matrix",sep="")

write.table(common.matrixA.cis,file=cisA.output.name,quote=FALSE,sep="\t")
write.table(common.matrixB.cis,file=cisB.output.name,quote=FALSE,sep="\t")
write.table(common.matrixA.trans,file=transA.output.name,quote=FALSE,sep="\t")
write.table(common.matrixB.trans,file=transB.output.name,quote=FALSE,sep="\t")
