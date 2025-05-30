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
#'         Should be fitted using \code{\link{estimateBestQ}} if unknown.
#'         If a named vector is provided it is used for the different
#'         splicing types.
#' @param implementation The method that should be used to correct for 
#' confounders. 
#' @param type The type of PSI (jaccard, psi5, psi3 or theta for theta/splicing 
#' efficiency)
#' @param iterations The maximal number of iterations. When the autoencoder has 
#' not yet converged after these number of iterations, the fit stops anyway.
#' @param BPPARAM A BiocParallel object to run the computation in parallel
#' @param correction Deprecated. The name changed to implementation. 
#' @param ... Additional parameters passed on to the internal fit function
#' 
#' @return FraserDataSet
#' @examples
#' # set default parallel backend
#' register(SerialParam())
#' 
#' # preprocessing
#' fds <- createTestFraserDataSet()
#' 
#' # filtering not expressed introns
#' fds <- calculatePSIValues(fds)
#' fds <- filterExpressionAndVariability(fds)
#'
#' # Run the full analysis pipeline: fits distribution and computes p values
#' fds <- FRASER(fds, q=2, implementation="PCA")
#'
#' # afterwards, the fitted fds-object can be saved and results can 
#' # be extracted and visualized, see ?saveFraserDataSet, ?results and 
#' # ?plotVolcano
#'  
#' # The functions run inside the FRASER function can also be directly 
#' # run themselves. 
#' # To directly run the fit function:
#' fds <- fit(fds, implementation="PCA", q=2, type="jaccard")
#' 
#' # To directly run the nomial and adjusted p value and z score 
#' # calculation, the following functions can be used:
#' fds <- calculatePvalues(fds, type="jaccard")
#' head(pVals(fds, type="jaccard"))
#' fds <- calculatePadjValues(fds, type="jaccard", method="BY")
#' head(padjVals(fds, type="jaccard"))
#' fds <- calculateZscore(fds, type="jaccard")
#' head(zScores(fds, type="jaccard")) 
#' 
#' # example of restricting FDR correction to subsets of genes of interest
#' genesOfInterest <- list("sample1"=c("TIMMDC1"), "sample2"=c("MCOLN1"))
#' fds <- calculatePadjValues(fds, type="jaccard", 
#'                  subsets=list("exampleSubset"=genesOfInterest))
#' padjVals(fds, type="jaccard", subsetName="exampleSubset")
#' padjVals(fds, type="jaccard", level="gene", subsetName="exampleSubset")
#' fds <- calculatePadjValues(fds, type="jaccard", 
#'                  subsets=list("anotherExampleSubset"=c("TIMMDC1")))
#' padjVals(fds, type="jaccard", subsetName="anotherExampleSubset")
#' 
#' # only adding FDR corrected pvalues on a subset without calculating 
#' # transcriptome-wide FDR again:
#' fds <- calculatePadjValuesOnSubset(fds, genesToTest=genesOfInterest, 
#'          subsetName="setOfInterest", type="jaccard")
#' padjVals(fds, type="jaccard", subsetName="setOfInterest")
#' 
#' @seealso \code{\link[FRASER]{fit}}
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @author Ines Scheller \email{scheller@@in.tum.de}
#' 
#' @rdname FRASER
#' @name FRASER
NULL

#' @describeIn FRASER This function runs the default FRASER pipeline combining 
#' the beta-binomial fit, the computation of Z scores and p values as well as 
#' the computation of delta-PSI values.
#' @export
FRASER <- function(fds, q, type=fitMetrics(fds), 
                   implementation=c("PCA", "PCA-BB-Decoder", "AE-weighted", 
                                    "AE", "BB"), 
                    iterations=15, BPPARAM=bpparam(), correction, 
                    subsets=NULL, ...){
    # Check input
    implementation <- match.arg(implementation)
    if (!missing("correction")){
        warning("The argument correction is deprecated. Use implementation ",
                "instead.")
        implementation <- correction
    }
    checkFraserDataSet(fds)

    # compute stats if needed
    if(isFALSE(checkCountData(fds, FALSE))){
        fds <- calculatePSIValues(fds)
    }

    # fit each splicing type separately
    for(i in type){

        # get type specific q
        if(missing(q)){
            currQ <- bestQ(fds, i)
        } else if(i %in% names(q)){
            currQ <- q[i]
        } else {
            warning("You provided one q for all metric. ",
                    "We recommend to estimate it per metric, ",
                    "which will give better and more reliable results.")
            currQ <- q
        }

        # fit autoencoder
        message("\n", date(), ": Fit step for: '", i, "'.")
        fds <- fit(object=fds, implementation=implementation, q=currQ,
                iterations=iterations, type=i, BPPARAM=BPPARAM, ...)

        message(date(), ": Compute p values for: '", i, "'.")
        fds <- calculatePvalues(fds, type=i)

        message(date(), ": Adjust p values for: '", i, "'.")
        fds <- calculatePadjValues(fds, type=i, subsets=subsets)

        # message(date(), ": Compute Z scores for: '", i, "'.")
        # fds <- calculateZscore(fds, type=i)
    }

    # return final analysis
    return(fds)
}
