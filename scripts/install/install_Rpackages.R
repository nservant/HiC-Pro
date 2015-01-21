## Nicolas Servant
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
