#'
#' Create example data sets for FraseR
#'
#' Creates an example data set from files with example counts. 
#' \code{makeFittedExampleFraseRDataSet} additionally runs the FRASER 
#' pipeline.
#' 
#' @param workingDir directory where to store HDF5 and RDS files. Defaults to
#'                the current tempory R session folder.
#' @param rerun Defaults to \code{FALSE}. If set to \code{TRUE} it reruns the
#'                full fit of the model.
#' 
#' @return FraseRDataSet
#'
#' @examples
#'   fds <- makeExampleFraseRDataSet()
#'   fds
#'   
#' @rdname makeExampleFraseRDataSet
#' @aliases makeExampleFraseRDataSet makeFittedExampleFraseRDataSet
#' @export
makeExampleFraseRDataSet <- function(workingDir=tempdir()){
    # example sample annoation
    sampleTable <- fread(system.file(
        "extdata", "sampleTable_countTable.tsv", package="FRASER", mustWork=TRUE))
    
    # get raw counts 
    junctionCts   <- fread(system.file("extdata", "raw_junction_counts.tsv.gz",
                                       package="FRASER", mustWork=TRUE))
    
    spliceSiteCts <- fread(system.file("extdata", "raw_site_counts.tsv.gz",
                                       package="FRASER", mustWork=TRUE))
    
    # create FraseR object
    fds <- FraseRDataSet(colData=sampleTable, junctions=junctionCts,
                         spliceSites=spliceSiteCts, workingDir=workingDir,
                         name="Example_Dataset")
    return(fds)
}

#' @rdname makeExampleFraseRDataSet
#' @export
makeFittedExampleFraseRDataSet <- function(workingDir=tempdir(), rerun=FALSE){
    # check if file exists already
    hdf5Files <- file.path(workingDir, "savedObjects", "Example_Dataset", 
                           "fds-object.RDS")
    if(all(file.exists(hdf5Files))){
        if(isFALSE(rerun)){
            fds <- loadFraseRDataSet(workingDir, name="Example_Dataset")
            if(all(paste0(c("zScores", "pajdBetaBinomial", "predictedMeans"),
                          "_", rep(psiTypes, 3)) %in% assayNames(fds))){
                message(date(), ": Use existing cache data.")
                return(fds)
            }
        }
        cleanCache(makeExampleFraseRDataSet(workingDir), all=TRUE)
    }
    
    # get test sample annotation
    fds <- makeExampleFraseRDataSet()
    
    # filter expression
    fds <- calculatePSIValues(fds)
    fds <- filterExpression(fds, minExpressionInOneSample=5, minDeltaPsi=0, 
                            quantileMinExpression=0)
    
    # run FraseR pipeline
    fds <- FraseR(fds, q=2, iterations=2)
    
    # annotate it
    fds <- annotateRanges(fds)
    
    # save data for later 
    dontWriteHDF5(fds) <- FALSE
    devNULL <- saveFraseRDataSet(fds)
    dontWriteHDF5(fds) <- TRUE
    
    # return a FraseRDataSet object
    return(fds)
}
