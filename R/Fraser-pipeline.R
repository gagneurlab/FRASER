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
#'   plotSampleResults(fds, 'sample1')
#'
FraseR <- function(fds, q, NcpuPerSample=1, ...){

    # Check input
    stopifnot(class(fds) == "FraseRDataSet")

    # count data
    if(!"rawCountsJ" %in% assayNames(fds))
        fds <- countRNAData(fds, NcpuPerSample=NcpuPerSample)

    # calculate PSI values
    fds <- calculatePSIValues(fds)


    # fit autoencoder
    if(missing(q)){
        warning("Please provide a fitted q to get better results!")
        q <- ceiling(ncol(fds)/10)
    }
    for(pt in psiTypes)
        fds <- fit(fds, q=q, type=pt, ...)

    # calculate ZScores
    for(pt in psiTypes)
        fds <- calculateZscore(fds, type=pt)

    # calculte P-values
    for(pt in psiTypes)
        fds <- calculatePvalues(fds, type=pt)

    # adjust pvalues
    for(pt in psiTypes)
        fds <- calculatePadjValues(fds, type=pt)

    # return final analysis
    return(fds)
}
