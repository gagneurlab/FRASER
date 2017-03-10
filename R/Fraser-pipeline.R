##
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains the standard pipeline for a FraseR analysis
##
#'
#' Find rare splicing events within RNA-seq data
#'
#' @description This function performs a default analysis of RNA-seq data
#'
#' @param settings A FraseRSetting object with all the information
#'             how and what to count
#' @param internBPPARAM A BiocParallel param object to configure the
#'             parallel backend of the internal loop for counting
#'
#' @return FraseRDataSet
#' @export
#' @examples
#'   fds <- FraseR(createTestFraseRSettings())
#'   plotSampleResults(fds)
FraseR <- function(settings, internBPPARAM=SerialParam()){

    # Check input
    stopifnot(class(settings) == "FraseRSettings")
    stopifnot(is(internBPPARAM, "BiocParallelParam"))

    # count data
    fds <- countRNAData(settings, internBPPARAM=internBPPARAM)

    # calculate PSI values
    fds <- calculatePSIValues(fds)

    # calculate ZScores
    fds <- calculateZScores(fds)

    # calculte P-values
    fds <- calculatePValues(fds)

    # return final analysis
    return(fds)
}
