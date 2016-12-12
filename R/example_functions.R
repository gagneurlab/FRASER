#'
#' Create a test case dataset to be used in the vignette and
#' to explore the functionallity of the FraseR package.
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
    
    # return a FraseRSettings object
    return(FraseRSettings(sampleData=sampleTable))
}