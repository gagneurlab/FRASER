
.getReadTypeFromPsiType <- function(psiType){
    stopifnot(length(psiType) == 1 && is(psiType, "character"))
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
