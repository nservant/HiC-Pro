## HiC-Pro                                                                                                           ## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Install missing R packages
##


##########################################
## CRAN

## RColorBrewer
if (!require(RColorBrewer)){
install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
}
##ggplot2
if (!require(ggplot2)){
install.packages("ggplot2", repos="http://cran.us.r-project.org")
}
##grid
if (!require(grid)){
install.packages("grid", repos="http://cran.us.r-project.org")
}
