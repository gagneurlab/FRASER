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
#' @return a FraseRDataSet object which contains a test case 
#' 
#' @examples
#' fds <- createTestFraseRSettings()
#' fds
#' 
#' fds <- createTestFraseRSettings("FRASER-data")
#' fds
#' 
#' fds <- createTestFraseRDataSet()
#' fds
#' 
#' @rdname createTestFraseRDataSet
#' @aliases createTestFraseRSettings createTestFraseRDataSet
#' @export
createTestFraseRSettings <- function(workingDir=tempdir()){

    # get sample data table
    sampleTable <- fread(system.file(
            "extdata", "sampleTable.tsv", package="FraseR", mustWork=TRUE))

    # convert it to a bamFile list
    bamFiles <- system.file(
        sampleTable[,bamFile], package="FraseR", mustWork=TRUE)
    sampleTable[,bamFile:=bamFiles]

    # TODO remove after example data update
    if("group" %in% colnames(sampleTable)){
        setnames(sampleTable, "group", "condition")
    }

    # check that NHDF is NA group
    sampleTable[gene=='NHDF', condition:=NA]
    
    # create FraseR object
    fds <- FraseRDataSet(colData=sampleTable, workingDir=workingDir)

    # dont use hdf5 for example data set
    dontWriteHDF5(fds) <- TRUE

    # return a FraseRSettings object
    fds
}


#' @rdname createTestFraseRDataSet
#' @export
createTestFraseRDataSet <- function(workingDir=tempdir(), rerun=FALSE){
    # check if file exists already
    padjHDF5 <- file.path(workingDir, "savedObjects", "Data_Analysis", 
            "pajdBetaBinomial_psiSite.h5")
    if(file.exists(padjHDF5)){
        if(isTRUE(rerun)){
            cleanCache(createTestFraseRSettings(workingDir), all=TRUE)
        } else {
            message(date(), ": Use existing cache data.")
            return(loadFraseRDataSet(workingDir, name="Data_Analysis"))
        }
    }
    
    # get test sample annotation
    fds <- createTestFraseRSettings(workingDir)
    
    # count data
    fds <- countRNAData(fds)
    
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

