#'
#' FRASER: Find RAre Splicing Events in RNA-seq data
#'
#' This help page describes the FRASER function which can be used run the 
#' default FRASER pipeline. This pipeline combines the beta-binomial fit, the 
#' computation of Z scores and p values as well as the computation of delta-PSI 
#' values.
#' 
#' All computed values are returned as an FraserDataSet object. To have
#' more control over each analysis step, one can call each function separately.
#' \itemize{
#'     \item \code{fit} to control for confounding effects and fit the beta 
#'     binomial model parameters
#'     \item \code{calculatePvalues} to calculate the nominal p values
#'     \item \code{calculatePadjValues} to calculate adjusted p values (per 
#'     sample)
#'     \item \code{calculateZscore} to calculate the Z scores
#' }
#' 
#' Available methods to correct for the confounders are currently: a denoising 
#' autoencoder with a BB loss ("AE" and "AE-weighted"), PCA ("PCA"), a hybrid 
#' approach where PCA is used to fit the latent space and then the decoder of 
#' the autoencoder is fit using the BB loss ("PCA-BB-Decoder"). Although not 
#' recommended, it is also possible to directly fit the BB distrbution to the 
#' raw counts ("BB"). 
#'
#' @inheritParams countRNA
#' @param q The encoding dimensions to be used during the fitting proceadure.
#'         Should be fitted using \code{\link{optimHyperParams}} if unknown.
#'         If a named vector is provided it is used for the different
#'         splicing types.
#' @param implementation The method that should be used to correct for 
#' confounders. 
#' @param type The type of PSI (psi5, psi3 or theta for theta/splicing 
#' efficiency)
#' @param iterations The maximal number of iterations. When the autoencoder has 
#' not yet converged after these number of iterations, the fit stops anyway.
#' @param BPPARAM A BiocParallel object to run the computation in parallel
#' @param correction Deprecated. The name changed to implementation. 
#' @param ... Additional parameters passed on to the internal fit function
### Additional parameters of the internal fit function:
#' @param rhoRange Defines the range of values that rho parameter from the 
#' beta-binomial distribution is allowed to take. For very small values of rho, 
#' the loss can be instable, so it is not recommended to allow rho < 1e-8. 
#' @param weighted If TRUE, the weighted implementation of the autoencoder is 
#' used
#' @param noiseAlpha Controls the amount of noise that is added for the 
#' denoising autoencoder.
#' @param convergence The fit is considered to have converged if the difference 
#' between the previous and the current loss is smaller than this threshold.
#' @param verbose Controls the level of information printed during the fit.
#' @param minDeltaPsi Minimal delta psi of an intron to be be considered a 
#' variable intron. 
#' @param initialize If FALSE and a fit has been previoulsy run, the values 
#' from the previous fit will be used as initial values. If TRUE, 
#' (re-)initialization will be done. 
#' @param control List of control parameters passed on to optim().
#' @param nSubset The size of the subset to be used in fitting if subsetting is
#' used.
#' 
#' @return FraserDataSet
#' @examples
#'    # preprocessing
#'    fds <- createTestFraserDataSet()
#'   
#'    ### when running FRASER on a real dataset, one should run the following 
#'    ### two commands first (not run here to make the example run faster):
#'    # fds <- calculatePSIValues(fds)
#'    # fds <- filterExpressionAndVariability(fds)
#'
#'    # Run the full analysis pipeline: fits distribution and computes p values
#'    fds <- FRASER(fds, q=2, implementation="PCA")
#'
#'    # afterwards, the fitted fds-object can be saved and results can 
#'    # be extracted and visualized, see ?saveFraserDataSet, ?results and 
#'    # ?plotVolcano
#'    
#'    ### The functions run inside the FRASER function can also be directly 
#'    ### run themselves. 
#'    ### To directly run the fit function:
#'    # fds <- fit(fds, implementation="PCA", q=2, type="psi5")
#'    
#'    ### To directly run the nomial and adjusted p value and z score 
#'    ### calculation, the following functions can be used:
#'    # fds <- calculatePvalues(fds, type="psi5")
#'    # head(pVals(fds, type="psi5"))
#'    # fds <- calculatePadjValues(fds, type="psi5", method="BY")
#'    # head(padjVals(fds, type="psi5"))
#'    # fds <- calculateZscore(fds, type="psi5")
#'    # head(zScores(fds, type="psi5")) 
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @rdname FRASER
#' @name FRASER
NULL

#' @describeIn FRASER This function runs the default FRASER pipeline combining 
#' the beta-binomial fit, the computation of Z scores and p values as well as 
#' the computation of delta-PSI values.
#' @export
FRASER <- function(fds, q, implementation=c("PCA", "PCA-BB-Decoder", 
                                        "AE-weighted", "AE", "BB"), 
                    iterations=15, BPPARAM=bpparam(), correction, ...){
    implementation <- match.arg(implementation)
    if (!missing("correction")){
        warning("The argument correction is deprecated. Use implementation ",
                "instead.")
        implementation <- correction
    }
    # Check input
    checkFraserDataSet(fds)

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
        fds <- fit(fds, implementation=implementation, q=currQ,
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
