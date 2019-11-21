#'
#' Create example data sets for FraseR
#'
#' Creates an example data set from a file or simulates a data set based
#' on  random split read counts following a beta-binomial distribution with
#' injected outliers.
#' 
#' @return FraseRDataSet
#'
#' @examples
#'   fds <- makeExampleFraseRDataSet()
#'   fds
#'
#' @export
makeExampleFraseRDataSet <- function(){
    dir <- system.file(package="FraseR")
    anno <- fread(file.path(dir, "extdata", "sampleTable.tsv"))
    anno[,bamFile:=file.path(dir,
            paste(strsplit(bamFile, "/")[[1]], collapse=.Platform$file.sep)),
            by=bamFile]

    fds <- FraseRDataSet(colData=anno)
    fds <- countRNAData(fds)
    fds <- calculatePSIValues(fds)
    return(fds)
}
