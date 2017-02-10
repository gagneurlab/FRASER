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
    sampleTable[,bamFile:=unlist(BamFileList(system.file(
                sampleTable[,bamFile], package="FraseR", mustWork=TRUE
    )))]
    
    # check that NHDF is NA group 
    sampleTable[gene=='NHDF',group:=NA]
    
    # return a FraseRSettings object
    return(FraseRSettings(
        sampleData=sampleTable, 
        outputFolder=file.path(Sys.getenv("HOME"), "FraseR")
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
    dataset <- calculateSitePSIValue(dataset)
    
    # return a FraseRDataSet object
    return(dataset)
}

