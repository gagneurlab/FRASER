#'
#' make an FraseR example data set
#' using an subset of the Kremer paper (200 junctions)
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
