# Reg3DIMS

Register sequential 2D IMS data loaded into `Cardinal` using `NiftyReg`.

`Reg3DIMS` extends [`Cardinal`](https://cardinalmsi.org)'s excellent data structures for 3D IMS experiments registered using [`RNiftyReg`](https://github.com/jonclayden/RNiftyReg). 

## Installation 

```R
#install devtools package is not already installed...
install.packages('devtools')

#install Reg3DIMS from github using devtools
library(devtools)
devtools::install_github('nhpatterson/Reg3DIMS')


#also install dependencies
install.packages('RNiftyReg')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBImage")
BiocManager::install("Cardinal")


```

## Markdown example
[3D pipeline on example Mouse Brain dataset](https://htmlpreview.github.io/?https://github.com/NHPatterson/Reg3DIMS/blob/master/markdown/Reg3DIMS.html)
