#'
#' Filtering FraseRDataSets
#'
#' @export
filterExpression <- function(fds, minExpressionInOneSample=20, quantile=0.05,
                    quantileMinExpression=1, minDeltaPsi=0.3, filter=TRUE,
                    BPPARAM=bpparam(), delayed=FALSE){

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
            rowQuantiles(as.matrix(ctsN5), probs=quantile) }
    f3 <- function(cts, ctsN3, quantile, ...) {
            rowQuantiles(as.matrix(ctsN3), probs=quantile) }
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

#'
#' Plot filter expression
#'
#' Histogram of the geometric mean per junction based on the filter status
#'
#' @export
plotFilterExpression <- function(fds, bins=200, alpha=0.75){
    cts    <- K(fds, "psi5")
    rowlgm <- exp(rowMeans(log(cts + 1)))
    dt <- data.table(
        value=rowlgm,
        passed=mcols(fds, type="j")[['passed']])
    colors <- brewer.pal(3, "Dark2")
    ggplot(dt, aes(value, fill=passed)) +
        geom_histogram(data=dt[passed==TRUE],
                aes(fill=colors[1]), alpha=alpha, bins=bins) +
        geom_histogram(data=dt[passed==FALSE],
                aes(fill=colors[2]), alpha=alpha, bins=bins) +
        scale_x_log10() +
        scale_y_log10() +
        scale_fill_manual(values=colors[1:2], name="Passed",
                labels=c("True", "False")) +
        xlab("Mean Junction Expression") +
        ylab("Count") +
        ggtitle("Expression filtering")
}
