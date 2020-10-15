# FRASER - Find RAre Splicing Events in RNA-seq

[![Build Status](https://travis-ci.com/c-mertes/FRASER.svg?branch=master)](https://travis-ci.com/c-mertes/FRASER)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/371x22cn6fipu7bw/branch/master?svg=true)](https://ci.appveyor.com/project/c-mertes/fraser/branch/master)
[![Version](https://img.shields.io/github/v/release/c-mertes/FRASER)](https://github.com/c-mertes/FRASER/releases)
[![Coverage status](https://codecov.io/gh/c-mertes/FRASER/branch/master/graph/badge.svg)](https://codecov.io/github/c-mertes/FRASER/branch/master)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/c-mertes/FRASER/blob/master/LICENSE)

FRASER is a tool to detect aberrant splicing events in RNA-seq data. The method is described in a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.18.866830v1) and available through [Bioconductor](http://bioconductor.org/packages/release/bioc/html/FRASER.html). It is also part of the [Detection of RNA Outlier Pipeline (DROP)](https://github.com/gagneurlab/drop). The DROP pipeline is described [here on Nature's protocol exchange](https://doi.org/10.21203/rs.2.19080/v1).
                                                                             
The FRASER framework and workflow aims to assist the diagnostics in the field of rare diseases where RNA-seq is performed to identify aberrant splicing defects. For a short tutorial on how to use FRASER on a dataset please use the [vignette](http://bioconductor.org/packages/release/bioc/vignettes/FRASER/inst/doc/FRASER.pdf) or our Colab tutorial at: [http://tinyurl.com/RNA-ASHG-colab](http://tinyurl.com/RNA-ASHG-colab). The Colab is based on a workshop that we presented at ASHG 2019/2020.

## Installation

`FRASER` is an R/Bioconductor software package requiring a running 
[R 3.6 version or higher](https://cran.r-project.org/).

The recommanded way of installing `FRASER` is to use `Bioconductor`:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages('BiocManager')
BiocManager::install('FRASER')
```

If you use an R version below `4.0.0` you have to install it from source. 
To this end, we will use `devtools` to install it. For this, you need a 
working development environment to compile the C++ code (see for 
details: [Windows](https://cran.r-project.org/bin/windows/Rtools/)
or [MacOS X](https://cran.r-project.org/bin/macosx/tools/)).

```
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages('devtools')

# latest development version
devtools::install_github('gagneurlab/FRASER', dependencies=TRUE)

# or a specific version of FRASER (search for tags on github)
devtools::install_github('gagneurlab/FRASER', ref='1.1.3', dependencies=TRUE)
```

If you have dependency issues while installing any package, please have a look
at the Troubleshooting section or submit an issue on 
[GitHub](https://github.com/gagneurlab/FRASER/issues).


### Toubleshooting

#### Missing libraries while compiling R packages

On some Linux distributions we need the developer libraries for compiling the R packages.

To install those packages, please run as administrator: 

For Ubuntu or debian based systems:
```
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libmysqld-dev
```

For centOS or RHEL based systems:
```
sudo yum install R-devel zlib-devel openssl-devel libcurl-devel libxml2-devel mariadb-devel
```
