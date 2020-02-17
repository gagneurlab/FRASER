#'
#' Filtering FraseRDataSets
#' 
#' This method can be used to filter out introns that are not reliably detected 
#' and to remove introns with no variablity between samples.
#' 
#' @inheritParams countRNA
#' @param minExpressionInOneSample The minimal read count in at least one 
#' sample that is required for an intron to pass the filter.
#' @param quantile Defines which quantile should be considered for the filter.
#' @param quantileMinExpression The minimum read count an intron needs to have 
#' at the specified quantile to pass the filter.
#' @param minDeltaPsi Only introns for which the maximal difference in the psi 
#' value of a sample to the mean psi of the intron is larger than this value 
#' pass the filter. 
#' @param filter If TRUE, a subsetted fds containing only the introns that 
#' passed all filters is returned. If FALSE, no subsetting is done and the 
#' information of whether an intron passed the filters is only stored in the 
#' mcols.
#' @param delayed If FALSE, count matrices will be loaded into memory, 
#' otherwise the function works on the delayedMatrix representations.
#'
#' @return FraseRDataSet
#' 
#' @examples
#' fds <- makeExampleFraseRDataSet()
#' fds <- calculatePSIValues(fds)
#' fds <- filterExpression(fds)
#'
#' @export
filterExpression <- function(fds, minExpressionInOneSample=20, quantile=0.05,
                    quantileMinExpression=1, minDeltaPsi=0, filter=TRUE,
                    BPPARAM=bpparam(), delayed=TRUE){

    # extract counts
    cts  <- K(fds, type="j")
    ctsN5 <- N(fds, type="psi5")
    ctsN3 <- N(fds, type="psi3")

    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        ctsN5 <- as.matrix(ctsN5)
        ctsN3 <- as.matrix(ctsN3)
    }

    # cutoff functions
    f1 <- function(cts, ...){
            rowMaxs(cts) }
    f2 <- function(cts, ctsN5, quantile, ...){
            rowQuantiles(ctsN5, probs=quantile) }
    f3 <- function(cts, ctsN3, quantile, ...) {
            rowQuantiles(ctsN3, probs=quantile) }
    f4 <- function(cts, ctsN3, ...) {
            psi <- cts/ctsN3
            rowMaxs(abs(psi - rowMeans2(psi, na.rm=TRUE)), na.rm=TRUE) }
    f5 <- function(cts, ctsN5, ...) {
            psi <- cts/ctsN5
            rowMaxs(abs(psi - rowMeans2(psi, na.rm=TRUE)), na.rm=TRUE) }

    funs <- c(maxCount=f1, quantileValue5=f2, quantileValue3=f3,
            maxDPsi3=f4, maxDPsi5=f5)

    # run it in parallel
    cutoffs <- bplapply(funs, function(f, ...) f(...), BPPARAM=BPPARAM,
            cts=cts, ctsN3=ctsN3, ctsN5=ctsN5, quantile=quantile)

    # add annotation to object
    for(n in names(cutoffs)){
        mcols(fds, type="j")[n] <- cutoffs[[n]]
    }
    mcols(fds, type="j")[['passed']] <-
            cutoffs$maxCount >= minExpressionInOneSample &
            (cutoffs$quantileValue5    >= quantileMinExpression |
                cutoffs$quantileValue3 >= quantileMinExpression) &
            (cutoffs$maxDPsi3     >= minDeltaPsi |
                cutoffs$maxDPsi5 >= minDeltaPsi)

    # filter if requested
    if(isTRUE(filter)){
        numFilt <- sum(mcols(fds, type="j")[['passed']])
        message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                ". This is ", signif(numFilt/length(fds)*100, 3),
                "% of the junctions"))
        fds <- fds[mcols(fds, type="j")[['passed']], by="psi5"]
    }

    validObject(fds)
    return(fds)
}
