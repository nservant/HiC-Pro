## Nicolas Servant
## Install R packages for HiC-Pro

p1<-require("RColorBrewer")
p2<-require("ggplot2")
p3<-require("grid")

if(p1 & p2 & p3){
     print("SUCCESS IN LOADING PACKAGES")
}else{
     stop("ERROR IN LOADING PACKAGES")	
}
	
