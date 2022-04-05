BTYPE <- "both"
NCPUS <- 3
Sys.setenv(MAKEFLAGS = "-j3")
BIOC_VERSION <- Sys.getenv("BIOC_VERSION")
START_TIME <- Sys.time()

print_log <- function(...){
    hash_line <- paste0(rep("#", 10), collapse="")
    message(paste0("\n", hash_line, "\n### ", date(), ": ", ..., "\n", hash_line, "\n"))
}

installIfReq <- function(p, type=BTYPE, Ncpus=NCPUS, ...){
    for(j in p){
        if(!requireNamespace(j, quietly=TRUE)){
            print_log("Install ", j)
            BiocManager::install(j, type=type, Ncpus=Ncpus, ...)
        }
    }
}

# install Bioconductor
if(!requireNamespace("BiocManager", quietly=TRUE)){
    print_log("Install BiocManager")
    install.packages("BiocManager", Ncpus=NCPUS)
}
BiocManager::install("BiocVersion", version=BIOC_VERSION)

# because of https://github.com/r-windows/rtools-installer/issues/3
if("windows" == .Platform$OS.type){
    print_log("Install XML on windows ...")
    BTYPE <- "win.binary"
    installIfReq(p=c("XML", "xml2", "RSQLite", "progress", "tibble", "AnnotationDbi", "BiocCheck", "rtracklayer"))
    
    print_log("Install source packages only for windows ...")
    BiocManager::install(c("GenomeInfoDbData", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"), 
            type="both", version=BIOC_VERSION)
}

# install needed packages
# add testthat to pre installation dependencies due to: https://github.com/r-lib/pkgload/issues/89
for(p in c("getopt", "XML", "xml2", "testthat", "devtools", "covr",
            "roxygen2", "BiocCheck", "R.utils", "GenomeInfoDbData",
            "rtracklayer", "hms")){
    installIfReq(p=p, type=BTYPE, Ncpus=NCPUS)
}

# install package with its dependencies with a timeout due to travis 50 min
R.utils::withTimeout(timeout=2400, {
    try({
        print_log("Update packages")
        BiocManager::install(ask=FALSE, type=BTYPE, Ncpus=NCPUS, version=BIOC_VERSION)
 
        print_log("Install dev package")
        devtools::install(".", dependencies=TRUE, type=BTYPE)
    })
})

# to get FRASER session info
try({ library(FRASER) })
print(BiocManager::valid())
