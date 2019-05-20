#'
#' Filtering FraseRDataSets
#'
#' @export
filterExpression <- function(fds, minExpressionInOneSample=10, quantile=0.5,
                    quantileMinExpression=1, minDeltaPsi=0.025, filter=FALSE,
                    BPPARAM=bpparam(), delayed=FALSE){

    # extract counts
    cts  <- K(fds, type="j")
    cts5 <- N(fds, type="psi5")
    cts3 <- N(fds, type="psi3")

    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        cts5 <- as.matrix(cts5)
        cts3 <- as.matrix(cts3)
    }

    # cutoff functions
    f1 <- function(cts, ...)                 rowMaxs(cts)
    f2 <- function(cts, cts5, quantile, ...){
            rowQuantiles(cts + cts5, probs=quantile) }
    f3 <- function(cts, cts3, quantile, ...) {
            rowQuantiles(cts + cts3, probs=quantile) }
    f4 <- function(cts, cts3, ...) {
            psi <- cts/(cts + cts3)
            rowMaxs(abs(t(t(psi) - rowMeans2(psi, na.rm=TRUE))), na.rm=TRUE) }
    f5 <- function(cts, cts5, ...) {
            psi <- cts/(cts + cts5)
            rowMaxs(abs(t(t(psi) - rowMeans2(psi, na.rm=TRUE))), na.rm=TRUE) }

    funs <- c(maxCount=f1, quantileValue5=f2, quantileValue3=f3,
            maxDPsi3=f4, maxDPsi5=f5)

    # run it in parallel
    cutoffs <- bplapply(funs, function(x, ...) x(...),
            cts=cts, cts3=cts3, cts5=cts5, quantile=quantile, BPPARAM=BPPARAM)

    # add annotation to object
    for(n in names(cutoffs)){
        mcols(fds, type="j")[n] <- cutoffs[[n]]
    }
    mcols(fds, type="j")['passed'] <-
            cutoffs$maxCount >= minExpressionInOneSample &
            (cutoffs$quantileValue5    >= quantileMinExpression |
                cutoffs$quantileValue3 >= quantileMinExpression) &
            (cutoffs$maxDPsi3     >= minDeltaPsi |
                 cutoffs$maxDPsi5 >= minDeltaPsi)

    # filter if requested
    if(filter==TRUE){
        numFilt <- sum(mcols(fds, type="j")[,'passed'])
        warning(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                ". This is ", signif(numFilt/length(fds)*100, 3),
                "% of the junctions"))
        fds <- fds[mcols(fds, type="j")[,'passed']]
    }

    validObject(fds)
    return(fds)
}
