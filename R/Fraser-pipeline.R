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
#' @param settings A FraseRDataSet object with all the information
#'             how and what to count
#' @param NcpuPerSample A BiocParallel param object to configure the
#'             parallel backend of the internal loop for counting
#'
#' @return FraseRDataSet
#' @export
#' @examples
#'   fds <- createTestFraseRSettings()
#'   fds <- FraseR(fds)
#'
#'   # save the final FraseR object
#'   saveFraseRDataSet(fds)
#'
#'   # finally visualize the results
#'   plotSampleResults(fds)
#'
FraseR <- function(settings=createTestFraseRSettings(), NcpuPerSample=1){

    # Check input
    stopifnot(class(settings) == "FraseRDataSet")
    # stopifnot(is(NcpuPerSample, "BiocParallelParam"))

    # count data
    fds <- countRNAData(settings, NcpuPerSample=NcpuPerSample)

    # calculate PSI values
    fds <- calculatePSIValues(fds)

    # calculate ZScores
    fds <- calculateZScores(fds)

    # calculte P-values
    fds <- calculatePValues(fds)

    # return final analysis
    return(fds)
}
