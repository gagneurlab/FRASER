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

getFraseR <- function(clean=FALSE, count=TRUE){
    fds <- NULL
    try({
        fds <- loadFraseRDataSet(getDir(), getName())
        if(clean == TRUE){
            cleanCache(fds, all=TRUE)
            fds <- getFraseR(count=count)
        }
    }, silent=TRUE)
    if(is.null(fds)) {
        fds <- createTestFraseRSettings()
        workingDir(fds) <- getDir()
        name(fds) <- getName()
        parallel(fds) <- MulticoreParam(4)

        if(count==TRUE){
            fds <- FraseR(settings=fds)
            fds <- saveFraseRDataSet(fds)
        }
    }
    return(fds)
}

invisible(getFraseR(TRUE, FALSE))
