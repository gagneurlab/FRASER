#'
#' Filtering FraseRDataSets
#'
#' @export
filterExpression <- function(fds, minExpressionInOneSample=5, quantile=0.7,
                    quantileMinExpression=1, filter=FALSE){

    # extract counts
    cts  <- as.matrix(counts(fds, type="j", side="of"))
    cts5 <- as.matrix(counts(fds, type="psi5", side="ot"))
    cts3 <- as.matrix(counts(fds, type="psi3", side="ot"))

    # cutoff functions
    f1 <- function() rowMaxs(cts)
    f2 <- function() rowQuantiles(cts + cts5, probs=quantile)
    f3 <- function() rowQuantiles(cts + cts3, probs=quantile)
    funs <- c(maxVal=f1, quantileVal5=f2, quantileVal3=f3)

    # run it in parallel
    cutoffs <- parallel::mclapply(funs, function(x) x(), mc.cores=3)

    # add annotation to object
    mcols(fds, type="j")['maxCount'] <- cutoffs$maxVal
    mcols(fds, type="j")['quantileValue5'] <- cutoffs$quantileVal5
    mcols(fds, type="j")['quantileValue3'] <- cutoffs$quantileVal3
    mcols(fds, type="j")['passed'] <-
            cutoffs$maxVal >= minExpressionInOneSample &
            (cutoffs$quantileVal5 >= quantileMinExpression |
                cutoffs$quantileVal3 >= quantileMinExpression)

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



filterExpressionDelayed <- function(fds, minExpressionInOneSample=5, quantile=0.7,
                                    quantileMinExpression=1, filter=FALSE){

    # extract counts
    cts  <- counts(fds, type="j", side="of")
    cts5 <- counts(fds, type="psi5", side="ot")
    cts3 <- counts(fds, type="psi3", side="ot")

    # cutoff functions
    # add annotation to object
    mcols(fds, type="j")['maxCount'] <- rowMaxs(cts)
    mcols(fds, type="j")['quantileValue5'] <-
        rowQuantiles(cts + cts5, probs=quantile, drop=FALSE)
    mcols(fds, type="j")['quantileValue3'] <-
        rowQuantiles(cts + cts3, probs=quantile, drop=FALSE)
    mcols(fds, type="j")['passed'] <-
        as.matrix(mcols(fds, type="j")['maxCount']) >= minExpressionInOneSample &
        (as.matrix(mcols(fds, type="j")['quantileValue5']) >= quantileMinExpression |
             as.matrix(mcols(fds, type="j")['quantileValue3']) >= quantileMinExpression)

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



