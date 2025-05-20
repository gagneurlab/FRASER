# FRASER - Find RAre Splicing Events in RNA-seq

[![R-CMD-check-bioc](https://github.com/gagneurlab/FRASER/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/gagneurlab/FRASER/actions/workflows/check-bioc.yml)
[![Version](https://img.shields.io/github/v/release/c-mertes/FRASER)](https://github.com/c-mertes/FRASER/releases)
[![Coverage status](https://codecov.io/gh/c-mertes/FRASER/branch/master/graph/badge.svg)](https://codecov.io/github/c-mertes/FRASER/branch/master)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/c-mertes/FRASER/blob/master/LICENSE)

FRASER is a tool to detect aberrant splicing events in RNA-seq data. The method is described and published in [Nature Communications](https://doi.org/doi:10.1038/s41467-020-20573-7) and available through [Bioconductor](https://doi.org/doi:10.18129/B9.bioc.FRASER). It is also part of the [Detection of RNA Outlier Pipeline (DROP)](https://github.com/gagneurlab/drop). DROP is described and published in [Nature Protocols](https://doi.org/doi:10.1038/s41596-020-00462-5).
                                                                             
The FRASER framework and workflow aims to assist the diagnostics in the field of rare diseases where RNA-seq is performed to identify aberrant splicing defects. For a short tutorial on how to use FRASER on a dataset please use the [vignette](http://bioconductor.org/packages/release/bioc/vignettes/FRASER/inst/doc/FRASER.pdf) or our Colab tutorial at: [http://tinyurl.com/RNA-ASHG-colab](http://tinyurl.com/RNA-ASHG-colab). The Colab is based on a workshop that we presented at ASHG 2019/2020.

Please cite our method paper FRASER if you use it in a publication:

> Mertes, C., Scheller, I.F., Yépez, V.A. *et al.* Detection of aberrant splicing events in RNA-seq data using FRASER. *Nat Commun* **12**, 529 (2021). https://doi.org/10.1038/s41467-020-20573-7

or if you use FRASER2:

> Scheller, I.F., Lutz, K., Mertes, C *et al.* Improved detection of aberrant splicing with FRASER 2.0 and the intron Jaccard index. *Am Jrnl Hum Genet* **110**, 12 (2023). https://doi.org/10.1016/j.ajhg.2023.10.014

## What's new

FRASER 2.0, an improved version of FRASER, is now available and used by default (version 1.99.0 and above).
FRASER 2.0 uses the Intron Jaccard Index as its splice metric instead of FRASER's 
previous three metrics along with some other parameter optimizations of pseudocounts, 
filtering settings and default delta cutoff. 
 
To change the splice metric, set `fitMetrics(fds)` to one or more of the metrics 
specified in `FRASER::psiTypes`. For FRASER 2.0 and the Intron Jaccard Index, the 
new default delta cutoff is 0.1 instead of the previous value of 0.3. When using 
the 3 previous metrics, the delta cutoff should be set manually to 0.3 
during results extraction, e.g. `results(fds, deltaPsiCutoff=0.3,...)`.

## Installation

`FRASER` is an R/Bioconductor software package requiring a running 
[R 3.6 version or higher](https://cran.r-project.org/).

The recommended way of installing `FRASER` is through `Bioconductor`:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages('BiocManager')
BiocManager::install('FRASER')
```

If you use an R version below `4.0.0`, you have to install it from the source. 
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


### Troubleshooting

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

On some Linux distributions, we need the developer libraries for compiling the R packages.

To install those packages, please run as administrator: 

For Ubuntu or Debian-based systems:
```
sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev zlib1g-dev libmysqld-dev
```

For centOS or RHEL based systems:
```
sudo yum install R-devel zlib-devel openssl-devel libcurl-devel libxml2-devel mariadb-devel
```
