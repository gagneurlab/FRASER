#'
#' Create a test dataset
#' 
#' Create a test case dataset based on the test sample annotation to be 
#' used in the vignette and to explore the functionallity of the 
#' FRASER package. Dependent on the request only the sample annotation 
#' or a full fitted model is returned.
#' 
#' @param workingDir directory where to store HDF5 and RDS files. Defaults to
#'                the current tempory R session folder.
#' @param rerun Defaults to \code{FALSE}. If set to \code{TRUE} it reruns the
#'                full fit of the model.
#' @return a FraserDataSet object which contains a test case 
#' 
#' @examples
#' fds <- createTestFraserSettings()
#' fds
#' 
#' fds <- createTestFraserDataSet()
#' fds
#' 
#' @rdname createTestFraserDataSet
#' @aliases createTestFraserSettings createTestFraserDataSet
#' @export
createTestFraserSettings <- function(workingDir=tempdir()){

    # get sample data table
    sampleTable <- fread(system.file(
            "extdata", "sampleTable.tsv", package="FRASER", mustWork=TRUE))

    # convert it to a bamFile list
    bamFiles <- system.file(
        sampleTable[,bamFile], package="FRASER", mustWork=TRUE)
    sampleTable[,bamFile:=bamFiles]

    # TODO remove after example data update
    if("group" %in% colnames(sampleTable)){
        setnames(sampleTable, "group", "condition")
    }

    # check that NHDF is NA group
    sampleTable[gene=='NHDF', condition:=NA]
    
    # create FRASER object
    fds <- FraserDataSet(colData=sampleTable, workingDir=workingDir)

    # dont use hdf5 for example data set
    dontWriteHDF5(fds) <- TRUE

    # return a FraserSettings object
    fds
}


#' @rdname createTestFraserDataSet
#' @export
createTestFraserDataSet <- function(workingDir=tempdir(), rerun=FALSE){
    # check if file exists already
    hdf5Files <- file.path(workingDir, "savedObjects", "Data_Analysis", 
            "fds-object.RDS")
    if(all(file.exists(hdf5Files))){
        if(isFALSE(rerun)){
            fds <- loadFraserDataSet(workingDir, name="Data_Analysis")
            if(all(paste0(c("zScores", "padjBetaBinomial", "predictedMeans"),
                        "_", rep(psiTypes, 3)) %in% assayNames(fds))){
                message(date(), ": Use existing cache data.")
                return(fds)
            }
        }
        cleanCache(createTestFraserSettings(workingDir), all=TRUE)
    }
    
    # get test sample annotation
    fds <- createTestFraserSettings(workingDir)
    
    # count data
    fds <- countRNAData(fds, filter=FALSE)
    
    # filter expression
    fds <- calculatePSIValues(fds)
    fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=5, 
            minDeltaPsi=0, quantileMinExpression=0)
    
    # run FRASER pipeline
    fds <- FRASER(fds, q=2, iterations=2)
    
    # annotate it
    suppressMessages({ fds <- annotateRangesWithTxDb(fds) })
    
    # save data for later 
    fds <- saveFraserDataSet(fds)
    
    # return a FraserDataSet object
    return(fds)
}

