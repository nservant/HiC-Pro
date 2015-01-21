## HiC-Pro                                                                                                           ## Copyleft 2015 Institut Curie                                                                                      
## Author(s): Nicolas Servant                                                                                        
## Contact: nicolas.servant@curie.fr                                                                                 
## This software is distributed without any guarantee under the terms of the GNU General                             
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Check R packages import
##                                                                                                                   

p1<-require("RColorBrewer")
p2<-require("ggplot2")
p3<-require("grid")

if(p1 & p2 & p3){
     print("SUCCESS IN LOADING PACKAGES")
}else{
     stop("ERROR IN LOADING PACKAGES")	
}
	
