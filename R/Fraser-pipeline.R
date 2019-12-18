#'
#' FraseR: Find RAre Splicing Events in RNA-seq data
#'
#' The FraseR function runs the default FraseR pipeline combinig the fit,
#' the computation of Z scores and p values as well as the delta-PSI values.
#' All computed values are returned as an FraseRDataSet object. To have
#' more control over each analysis step, one can call each function separately.
#'
#' * controlForConfounders to control for confounding effects
#' * fitParams to fit the additional beta binomial model parameters
#'         (only needed if the autoencoder is not used)
#' * computePvalues to calculate the nominal and adjusted p values
#' * computeZscores to calculate the Z scores
#' * computeDeltaPsi to calculate the delta PSI values
#'
#' @inheritParams fit
#' @param fds A FraseRDataSet object
#' @param q The encoding dimensions to be used during the fitting proceadure.
#'         If a named vector is provided it is used for the different
#'         splicing types.
#' @param BPPARAM A BiocParallel object to run the computation in parallel
#' @param ... Additional parameters passed on to the internal fit function
#'
#' @return FraseRDataSet
#' @examples
#'   # preprocessing
#'   fds <- makeExampleFraseRDataSet()
#'   fds <- calculatePSIValues(fds)
#'   fds <- filterExpression(fds)
#'
#'   # Run analysis pipeline
#'   fds <- FraseR(fds, q=2)
#'   fds
#'
#'   # save the final FraseR object
#'   saveFraseRDataSet(fds)
#'
#'   # finally visualize the results
#'   plotVolcano(fds, 'sample1', 'psi5')
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
FraseR <- function(fds, q, correction="PCA", iterations=15,
                    BPPARAM=bpparam(), ...){

    # Check input
    checkFraseRDataSet(fds)

    # compute stats if needed
    if(isFALSE(checkCountData(fds, FALSE))){
        fds <- calculatePSIValues(fds)
    }

    # fit autoencoder
    if(missing(q)){
        warning("Please provide a fitted q to get better results!")
        q <- ceiling(ncol(fds)/10)
    }

    # fit each splicing type separately
    for(i in psiTypes){

        # get type specific q
        currQ <- q
        if(i %in% names(q)){
            currQ <- q[i]
        }

        message("\n", date(), ": Fit step for: '", i, "'.")
        fds <- fit(fds, correction=correction, q=currQ,
                iterations=iterations, type=i, BPPARAM=BPPARAM, ...)

        message(date(), ": Compute p values for: '", i, "'.")
        fds <- calculatePvalues(fds, type=i)

        message(date(), ": Adjust p values for: '", i, "'.")
        fds <- calculatePadjValues(fds, type=i)

        message(date(), ": Compute Z scores for: '", i, "'.")
        fds <- calculateZscore(fds, type=i)
    }

    # return final analysis
    return(fds)
}
