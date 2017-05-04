
#'
#' Create a test case dataset (sample information only)
#' to be used in the vignette and to explore the
#' functionallity of the FraseR package.
#'
#' @return a FraseRSettings object which contains a test case
#' @export
#' @examples
#'     createTestFraseRSettings()
#'
createTestFraseRSettings <- function(){
    # get sample data table
    sampleTable <- fread(system.file(
                "extdata", "sampleTable.tsv", package="FraseR", mustWork=TRUE
    ))

    # convert it to a bamFile list
    bamFiles <- system.file(sampleTable[,bamFile],
            package="FraseR", mustWork=TRUE)
    sampleTable[,bamFile:=bamFiles]

    # TODO remove after example data update
    if("group" %in% colnames(sampleTable)){
        setnames(sampleTable, "group", "condition")
    }

    # check that NHDF is NA group
    sampleTable[gene=='NHDF',condition:=NA]

    # return a FraseRSettings object
    return(FraseRDataSet(
        colData=sampleTable,
        workingDir=file.path(Sys.getenv("HOME"), "FraseR")
    ))
}

#'
#' Create a test case dataset based on the test sample annotation
#' filled with counts to be used in the vignette and
#' to explore the functionallity of the FraseR package.
#'
#' @return a FraseRDataSet object which contains a test case
#' @export
#' @examples
#'     createTestFraseRDataSet()
#'
createTestFraseRDataSet <- function(BPPARAM=NULL){

    # get test sample annotation
    settings <- createTestFraseRSettings()

    # set paralelle settings if given
    if(!is.null(BPPARAM)){
        settings@parallel <- BPPARAM
    }

    # count data
    dataset <- countRNAData(settings)

    # calculate PSI values
    dataset <- calculatePSIValues(dataset)
    dataset <- calculateZScores(dataset)

    # calculate pvalues
    dataset <- saveFraseRDataSet(dataset)
    dataset <- calculatePValues(dataset)

    # annotate it
    dataset <- annotateRanges(dataset)

    # return a FraseRDataSet object
    return(dataset)
}

