#'
#' FraseR: Find RAre Splicing Events in RNA-seq data
#'
#' The FraseR function runs the default FraseR pipeline combining the fit,
#' the computation of Z scores and p values as well as the delta-PSI values.
#' 
#' All computed values are returned as an FraseRDataSet object. To have
#' more control over each analysis step, one can call each function separately.
#' \itemize{
#'     \item \code{fit} to control for confounding effects and fit the beta 
#'     binomial model parameters, see \code{?fit} for details
#'     \item \code{calculatePvalues} to calculate the nominal p values, see 
#'     \code{?calculatePvalues} for details
#'     \item \code{calculatePadjValues} to calculate adjusted p values, see 
#'     \code{?calculatePadjValues} for details
#'     \item \code{calculateZscore} to calculate the Z scores, see 
#'     \code{?calculateZscore} for details
#' }
#' 
#' Available methods to correct for the confounders are currently: a denoising 
#' autoencoder with a BB loss ("AE" and "AE-weighted"), PCA ("PCA"), a hybrid 
#' approach where PCA is used to fit the latent space and then the decoder of 
#' the autoencoder is fit using the BB loss ("PCA-BB-Decoder"). Although not 
#' recommended, it is also possible to directly fit the BB distrbution to the 
#' raw counts ("BB"). 
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
#'    # preprocessing
#'    fds <- createTestFraseRDataSet()
#'   
#'    # when running FRASER on a real dataset, one should run the following 
#'    # two commands first (not run here to make the example run faster)
#'    # fds <- calculatePSIValues(fds)
#'    # fds <- filterExpressionAndVariability(fds)
#'
#'    # Run analysis pipeline: fits distribution and calculates p values
#'    fds <- FraseR(fds, q=2, correction="PCA")
#'    fds
#'
#'    # afterwards, the fitted fds-object can be saved and results can 
#'    # be extracted and visualized, see ?saveFraseRDataSet, ?results and 
#'    # ?plotVolcano
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
FraseR <- function(fds, q, correction=c("PCA", "PCA-BB-Decoder", "AE-weighted", 
                                        "AE", "BB"), iterations=15,
                    BPPARAM=bpparam(), ...){
    correction <- match.arg(correction)
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
