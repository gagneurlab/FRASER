#
# helper scripts for the testing step
#

getDir <- function(useHome=FALSE){
    if(useHome){
        file.path(Sys.getenv("HOME"),  "FraseR_test_that")
    } else {
        file.path(tempdir(), "FraseR_test_that")
    }
}

getName <- function(){
    "test_that"
}

getFraseRTestFile <- function(useHome=FALSE){
    file.path(getDir(useHome), getName(), "fds-object.RDS")
}

getFraseR <- function(clean=FALSE, useHome=FALSE){
    fds <- NULL
    try({
        fds <- readRDS(getFraseRTestFile(useHome=useHome))
        if(clean == TRUE){
            unlink(getFraseRTestFile(useHome=useHome))
            fds <- getFraseR()
        }
    }, silent=TRUE)
    if(is.null(fds)) {
        fds <- createTestFraseRSettings()
        workingDir(fds) <- getDir()
        name(fds) <- getName()
        verbose(fds) <- 0
        dontWriteHDF5(fds) <- TRUE

        if(.Platform$OS.type != "windows"){
            parallel(fds) <- MulticoreParam(4)
            register(parallel(fds))
        }

        fds <- countRNAData(fds)
        fds <- calculatePSIValues(fds)
        fds <- filterExpression(fds)
        fds <- FraseR(fds, q=2, correction="PCA")
        dir.create(dirname(getFraseRTestFile(useHome=useHome)), recursive=TRUE)
        saveRDS(fds, getFraseRTestFile(useHome=useHome))
    }
    return(fds)
}
