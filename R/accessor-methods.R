#'
#' @author Christian Mertes
#'
#' This file contains accessor functions for the FraseRDataSet object
#'

#'
#' the three different splice types
#'   * psi3: 3' splice site
#'   * psi5: 5' splice site
#'   * psiS: non splice reads vs spliced reads at splice site
#'
getPsiTypes <- function(){
    c("psi3", "psi5", "psiS")
}

getReadType <- function(psiType){
    stopifnot(isScalarCharacter(psiType))
    switch(psiType,
           psi3 = "splitReads",
           psi5 = "splitReads",
           psiS = "nonSplicedReads",
           stop("Did not recognice the given psi type: ", psiType)
    )
}

getPvalName <- function(){

}

fraserNames <- data.table(
    readType = c("splitReads", "splitReads", "nonSplicedReads"),
    psiType = c("psi3", "psi5", "sitePSI"),
    pvalName = c("pvalue_psi3", "pvalue_psi5", "pvalue_sitePSI")
)


fdsmcols <- function(fds, psiType){
    assays(fds)
}
