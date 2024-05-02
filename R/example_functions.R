#'
#' Create a test dataset
#' 
#' Create a test case dataset based on the test sample annotation to be 
#' used in the vignette and to explore the functionality of the 
#' FRASER package. Dependent on the request only the sample annotation 
#' or a full fitted model is returned.
#' 
#' @param workingDir Directory where to store HDF5 and RDS files. Defaults to
#'                \code{FRASER_output} in the current working directory.
#' @param rerun Defaults to \code{FALSE}. If set to \code{TRUE} it reruns the
#'                full fit of the model.
#' @param metrics The splice metrics that should be included in the test fds. 
#'              One or several of 'jaccard', 'psi5', 'psi3' or 'theta'.
#' @return A FraserDataSet object that contains a test case 
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
createTestFraserSettings <- function(workingDir="FRASER_output"){

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
    
    # set strand specificity column in colData
    sampleTable[, strand := 0L]
    
    # create FRASER object
    fds <- FraserDataSet(colData=sampleTable, workingDir=workingDir)

    # dont use hdf5 for example data set
    dontWriteHDF5(fds) <- TRUE

    # return a FraserSettings object
    fds
}


#' @rdname createTestFraserDataSet
#' @export
createTestFraserDataSet <- function(workingDir="FRASER_output", rerun=FALSE,
                                    metrics="jaccard"){
    metrics <- match.arg(metrics, several.ok=TRUE, 
                            choices=c("jaccard", "psi5", "psi3", "theta"))
    # check if file exists already
    hdf5Files <- file.path(workingDir, "savedObjects", "Data_Analysis", 
            "fds-object.RDS")
    if(all(file.exists(hdf5Files))){
        if(isFALSE(rerun)){
            fds <- loadFraserDataSet(workingDir, name="Data_Analysis")
            if(all(paste0(c("padjBetaBinomial", "predictedMeans"),
                        "_", rep(psiTypes, 3)) %in% assayNames(fds))){
                message(date(), ": Use existing cache data.")
                return(fds)
            }
        }
    }
    
    # get test sample annotation
    fds <- createTestFraserSettings(workingDir)
    
    # count data
    fds <- countRNAData(fds, filter=FALSE, recount=rerun)
    
    # filter expression
    fds <- calculatePSIValues(fds)
    fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=5, 
            minDeltaPsi=0, quantileMinExpression=0)
    
    # annotate it
    suppressMessages({ fds <- annotateRangesWithTxDb(fds) })
    
    # run FRASER pipeline
    fitMetrics(fds) <- metrics
    fds <- FRASER(fds, q=c(jaccard=2, psi5=2, psi3=2, theta=2), iterations=2)
    
    # save data for later 
    fds <- saveFraserDataSet(fds)
    
    # return a FraserDataSet object
    return(fds)
}

