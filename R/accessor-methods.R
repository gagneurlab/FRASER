#'
#' @author Christian Mertes
#'
#' This file contains accessor functions for the FraseRDataSet object
#'

#'
#' internal accessor for the non spliced reads object
#' @noRd
nonSplicedReads <- function(fds){
    stopifnot(class(fds) == "FraseRDataSet")
    return(fds@nonSplicedReads)
}

.getReadTypeFromPsiType <- function(psiType){
    stopifnot(length(psiType) == 1 && class(psiType) == "")
    readType <- switch(psiType,
           psi3="spliceSite",
           psi5="spliceSite",
           sitePSI="nonSplicedReads",
           stop("The given PSI type '", psiType, "' is not correct!")
    )
    return(readType)
}

fdsmcols <- function(fds, psiType){
    assays(fds)
}
