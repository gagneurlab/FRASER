#'
#' Create a test case dataset (sample information only)
#' to be used in the vignette and to explore the
#' functionallity of the FraseR package.
#'
#' @return a FraseRSettings object which contains a test case
#' @examples
#' fds <- createTestFraseRSettings()
#' fds
#'
#' @export
createTestFraseRSettings <- function(){

    # get sample data table
    sampleTable <- fread(system.file(
            "extdata", "sampleTable.tsv", package="FraseR", mustWork=TRUE
    ))

    # convert it to a bamFile list
    bamFiles <- system.file(
        sampleTable[,bamFile], package="FraseR", mustWork=TRUE
    )
    sampleTable[,bamFile:=bamFiles]

    # TODO remove after example data update
    if("group" %in% colnames(sampleTable)){
        setnames(sampleTable, "group", "condition")
    }

    # check that NHDF is NA group
    sampleTable[gene=='NHDF', condition:=NA]

    # create FraseR object
    fds <- FraseRDataSet(colData=sampleTable, parallel=bpparam(),
            workingDir=file.path(Sys.getenv("HOME"), "FraseR"))

    # dont use hdf5 for example data set
    dontWriteHDF5(fds) <- TRUE

    # return a FraseRSettings object
    fds
}

#'
#' Create a test case dataset based on the test sample annotation
#' filled with counts to be used in the vignette and
#' to explore the functionallity of the FraseR package.
#'
#' @return a FraseRDataSet object which contains a test case
#' @examples
#' fds <-  createTestFraseRDataSet()
#' fds
#'
#' @export
#'
createTestFraseRDataSet <- function(BPPARAM=bpparam()){

    # get test sample annotation
    fds <- createTestFraseRSettings()

    # set paralelle settings if given
    parallel(fds) <- BPPARAM

    # count data
    fds <- countRNAData(fds)

    # run FraseR pipeline
    fds <- FraseR(fds, q=2, iterations=2)

    # annotate it
    fds <- annotateRanges(fds)

    # return a FraseRDataSet object
    return(fds)
}

