##
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for calculating the statistics
## First it starts with calculating the Z-score for each
## site and then the p-values are calculated dependend on the
## given method in the setting file
##

#'
#' Calculate the zscore for each PSI value.
#'
#' @noRd
#' @return FraseRDataSet
#' @examples
#'   fds <- countRNAData(createTestFraseRSettings())
#'   fds <- calculatePSIValues(fds)
#'   fds <- calculateZScores(fds)
calculateZScores <- function(fds, type=psiTypes){

    # check input
    stopifnot(is(fds, "FraseRDataSet"))

    # calculate zscore for each psi type
    for(pt in type){
        fds <- calculateZScorePerDataSet(fds, pt)
    }

    return(fds)
}

#'
#' calculates the zscore for a given data type and a given psi type
#' and adds it directly to the dataset itself
#'
#' @noRd
calculateZScorePerDataSet <- function(fds, psiType){

    message(date(), ": Calculate the Zscore for ", psiType, " values ...")

    # get raw data and replace NA's with zeros
    psiVal <- assay(fds, psiType)

    # z = ( x - mean ) / sd
    rowmean <- rowMeans(psiVal, na.rm = TRUE)
    rowsd   <- apply(psiVal, 1, sd, na.rm = TRUE)
    zscores <- (psiVal - rowmean) / rowsd

    # use as.matrix to rewrite it as a new hdf5 array
    zScores(fds, type=psiType) <- as.matrix(zscores)

    return(fds)
}
