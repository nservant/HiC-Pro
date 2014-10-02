rm(list=ls())

### were additional arguments to R CMD BATCH given ?
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

matrix.import <- function(filename="",type="cis"){
  if(!file.exists(filename)){
    error.msg <- paste("can't fild matrix file -- ",filename ,sep="")
    stop(error.msg)
  }
  cur.matrix <- try(read.table(filename))
  if(class(cur.matrix)=='try-error'){
    if(type=="cis"){
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

cis.map.filename1 <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genome1,".matrix",sep="")
cis.map.filename2 <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome2,".",chrmstr2,"_",genome2,".matrix",sep="")
trans.map.filename <- paste(docDir,"/",sampleName,".",chrmstr1,"_",genome1,".",chrmstr2,"_",genome2,".matrix",sep="")

cis.matrix1 <-matrix.import(filename=cis.map.filename1,type="cis")
cis.matrix2 <-matrix.import(filename=cis.map.filename2,type="cis")
trans.matrix <-matrix.import(filename=trans.map.filename,type="trans")

cis.matrix <- cis.matrix1 + cis.matrix2
prob.matrix <- 1-trans.matrix/(cis.matrix+trans.matrix)

prob.matrix[is.infinite(prob.matrix)] <- 1
prob.matrix[is.na(prob.matrix)] <- 1

output.name <- paste(probDir,"/",sampleName,".",chrmstr1,".",chrmstr2,".cis.prob.matrix",sep="")
write.table(prob.matrix,file=output.name,quote=FALSE,sep="\t")



