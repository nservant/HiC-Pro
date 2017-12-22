## Nicolas Servant
## Install R packages for HiC-Pro

p1 <- require("RColorBrewer")
p2 <- require("ggplot2")
p2.1 <- packageVersion("ggplot2") >= '2.2.1'
p3<-require("grid")

if (! p2.1){
    stop("Error : ggplot2 version must be > 2.2.1 ")
}

if(p1 & p2 & p3){
     message("R packages loading with success")
}else{
     stop("Error in loading R packages - Please, be sure that the RColorBrewer, the ggplot2 and the grid packages are installed !")	
}
	
