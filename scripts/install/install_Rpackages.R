## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## Install R packages for HiC-Pro


##########################################
## CRAN

## RColorBrewer
if (!require(RColorBrewer)){
install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
}
if (!require(ggplot2)){
install.packages("ggplot2", repos="http://cran.us.r-project.org")
}
if (!require(grid)){
install.packages("grid", repos="http://cran.us.r-project.org")
}
