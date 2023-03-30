# FRASER - Find RAre Splicing Events in RNA-seq

[![Build](https://github.com/c-mertes/FRASER/workflows/Build/badge.svg?branch=master)](https://github.com/c-mertes/FRASER/actions?query=workflow%3ABuild)
[![Version](https://img.shields.io/github/v/release/c-mertes/FRASER)](https://github.com/c-mertes/FRASER/releases)
[![Coverage status](https://codecov.io/gh/c-mertes/FRASER/branch/master/graph/badge.svg)](https://codecov.io/github/c-mertes/FRASER/branch/master)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/c-mertes/FRASER/blob/master/LICENSE)

FRASER is a tool to detect aberrant splicing events in RNA-seq data. The method is described and published in [Nature Communications](https://doi.org/doi:10.1038/s41467-020-20573-7) and available through [Bioconductor](https://doi.org/doi:10.18129/B9.bioc.FRASER). It is also part of the [Detection of RNA Outlier Pipeline (DROP)](https://github.com/gagneurlab/drop). The DROP pipeline is described and published in [Nature Protocols](https://doi.org/doi:10.1038/s41596-020-00462-5).
                                                                             
The FRASER framework and workflow aims to assist the diagnostics in the field of rare diseases where RNA-seq is performed to identify aberrant splicing defects. For a short tutorial on how to use FRASER on a dataset please use the [vignette](http://bioconductor.org/packages/release/bioc/vignettes/FRASER/inst/doc/FRASER.pdf) or our Colab tutorial at: [http://tinyurl.com/RNA-ASHG-colab](http://tinyurl.com/RNA-ASHG-colab). The Colab is based on a workshop that we presented at ASHG 2019/2020.

Please cite our method paper if you use it in a publication:

> Mertes, C., Scheller, I.F., Yépez, V.A. *et al.* Detection of aberrant splicing events in RNA-seq data using FRASER. *Nat Commun* **12**, 529 (2021). https://doi.org/10.1038/s41467-020-20573-7

## What's new

FRASER 2.0, an improved version of FRASER that uses the Intron Jaccard Index as 
its splice metric instead of FRASER's previous three metrics, is now available 
and used by default (version 1.99.0 and above). The manuscript describing these 
changes in more detail will be available soon. 

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

#### Error in UseMethod("filter_")
When using FRASER with R3.6 one might observe the following error:

```
fds <- annotateRanges(fds)
# Error in UseMethod("filter_") :
#   no applicable method for 'filter_' applied to an object of class "c('tbl_SQLiteConnection', 'tbl_dbi', 'tbl_sql', 'tbl_lazy', 'tbl')"
```

To overcome this error one needs to upgrade BiocFileCache.

```
BiocManager::install("Bioconductor/BiocFileCache", ask=FALSE, update=FALSE)
```

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
