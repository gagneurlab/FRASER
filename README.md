# FRASER - Find RAre Splicing Events in RNA-seq
FRASER is a tool to detect aberrant splicing events in RNA-seq data. The method is currently avaiable as preprint [here](https://tinyurl.com/FRASER-paper).

[![Build Status](https://travis-ci.com/gagneurlab/FRASER.svg?branch=master)](https://travis-ci.com/gagneurlab/FRASER)
                                                                             
The FRASER framework and workflow aims to assist the diagnostics in the field of rare diseases where RNA-seq is performed to identify aberrant splicing defects. For a short tutorial on how to use FRASER on a dataset please use our Colab tutorial at: [http://tinyurl.com/RNA-ASHG-colab](http://tinyurl.com/RNA-ASHG-colab). This is based on a workshop we presented at ASHG 2019. 

## Installation

`FRASER` is an R software package requiring a running [R 3.4 version or higher](https://cran.r-project.org/).

We will use `devtools` to install it. For this, you need a working development environment to compile the C++ code (see for details: [Windows](https://cran.r-project.org/bin/windows/Rtools/)
or [MacOS X](https://cran.r-project.org/bin/macosx/tools/)).

```
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages('devtools')

# latest development version
devtools::install_github('gagneurlab/FRASER', dependencies=TRUE)
```

If you have dependency issues while installing any package, please have a look
at the Troubleshooting section or submit an issue on [GitHub](https://github.com/gagneurlab/FRASER/issues).


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
