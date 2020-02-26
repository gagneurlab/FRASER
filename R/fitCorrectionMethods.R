#'
#' Fit a beta-binomial distribution and correction for confounders
#' 
#' This method correct for confounders in the data and fits a beta-binomial 
#' distribution to the introns.
#' 
#' Available methods to correct for the confounders are currently: a denoising 
#' autoencoder with a BB loss ("AE" and "AE-weighted"), PCA ("PCA"), a hybrid 
#' approach where PCA is used to fit the latent space and then the decoder of 
#' the autoencoder is fit using the BB loss ("PCA-BB-Decoder"). Although not 
#' recommended, it is also possible to directly fit the BB distrbution to the 
#' raw counts ("BB"). 
#'
#' @inheritParams countRNA
#' @param correction The method that should be used to correct for confounders. 
#' @param q The dimension of the latent space. Should be fitted using 
#' \code{\link{optimHyperParams}} if unknown.
#' @param type The type of PSI (psi5, psi3 or psiSite for theta/splicing 
#' efficiency)
#' @param rhoRange Defines the range of values that rho parameter from the 
#' beta-binomial distribution is allowed to take. For very small values of rho, 
#' the loss can be instable, so it is not recommended to allow rho < 1e-8. 
#' @param weighted If TRUE, the weighted implementation of the autoencoder is 
#' used
#' @param noiseAlpha Controls the amount of noise that is added for the 
#' denoising autoencoder.
#' @param convergence The fit is considered to have converged if the difference 
#' between the previous and the current loss is smaller than this threshold.
#' @param iterations The maximal number of iterations. When the autoencoder has 
#' not yet converged after these number of iterations, the fit stops anyway.
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
#' @return FraseRDataSet
#' 
#' @examples 
#'   # generate toy data
#'   fds <- createTestFraseRDataSet()
#'   
#'   # fit
#'   fds <- fit(fds, correction="PCA", q=3, type="psi5")
#'   
#' @export
fit <- function(fds, correction=c("PCA", "PCA-BB-Decoder", "AE", "AE-weighted", 
                                "PCA-BB-full", "fullAE", "PCA-regression", 
                                "PCA-reg-full", "PCA-BB-Decoder-no-weights", 
                                "BB"),
                q, type="psi3", rhoRange=c(1e-8, 1-1e-8), 
                weighted=FALSE, noiseAlpha=1, convergence=1e-5, 
                iterations=15, initialize=TRUE, control=list(), 
                BPPARAM=bpparam(), nSubset=15000, verbose=FALSE, 
                minDeltaPsi=0.1){
    method <- match.arg(correction)

    # make sure its only in-memory data for k and n
    currentType(fds) <- type
    counts(fds, type=type, side="other", HDF5=FALSE) <- as.matrix(
            counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(
            counts(fds, type=type, side="ofInterest"))
    
    # check q is set
    if(correction != "BB" & (missing(q) | is.null(q))){
        stop("Please provide a q to define the size of the latent space!")
    }
    
    message(date(), ": Running fit with correction method: ", correction)
    fds <- switch(
        method,
        "AE"      = fitFraserAE(
            fds = fds, 
            q = q,
            type = type,
            noiseAlpha = noiseAlpha,
            rhoRange = rhoRange,
            lambda = 0,
            convergence = convergence,
            iterations = iterations,
            initialize = initialize,
            weighted = weighted,
            control = control,
            BPPARAM = BPPARAM,
            verbose = verbose,
            subset = TRUE,
            nrDecoderBatches = 1
        ),
        "AE-weighted" = fitFraserAE(
            fds = fds,
            q = q,
            type = type,
            noiseAlpha = noiseAlpha,
            nSubset = nSubset,
            rhoRange = rhoRange,
            lambda = 0,
            convergence = convergence,
            iterations = iterations,
            initialize = initialize,
            weighted = TRUE,
            control = control,
            BPPARAM = BPPARAM,
            verbose = verbose,
            subset = TRUE,
            nrDecoderBatches = 1
        ),
        "PCA-BB-Decoder" = fitFraserAE(
            fds = fds,
            q = q,
            type = type,
            noiseAlpha = noiseAlpha,
            nSubset = nSubset,
            rhoRange = rhoRange,
            lambda = 0,
            convergence = convergence,
            iterations = iterations,
            initialize = initialize,
            weighted = TRUE,
            control = control,
            BPPARAM = BPPARAM,
            verbose = verbose,
            subset = TRUE,
            nrDecoderBatches = 1,
            latentSpace = 'PCA'
        ),
        "PCA-BB-Decoder-no-weights" = fitFraserAE(
            fds = fds,
            q = q,
            type = type,
            noiseAlpha = noiseAlpha,
            nSubset = nSubset,
            rhoRange = rhoRange,
            lambda = 0,
            convergence = convergence,
            latentSpace = 'PCA',
            iterations = iterations,
            initialize = initialize,
            weighted = FALSE,
            control = control,
            BPPARAM = BPPARAM,
            verbose = verbose,
            subset = TRUE,
            nrDecoderBatches = 1
        ),
        "PCA-BB-full"       = fitFraserAE(
            fds = fds,
            q = q,
            type = type,
            noiseAlpha = noiseAlpha,
            rhoRange = rhoRange,
            lambda = 0,
            convergence = convergence,
            iterations = iterations,
            initialize = initialize,
            weighted = TRUE,
            control = control,
            BPPARAM = BPPARAM,
            verbose = verbose,
            subset = FALSE,
            nrDecoderBatches = 1,
            latentSpace = 'PCA'
        ),
        "PCA-reg-full"      = fitPCA(
            fds = fds,
            q = q,
            psiType = type,
            noiseAlpha = noiseAlpha,
            rhoRange = rhoRange,
            subset = FALSE,
            minDeltaPsi = minDeltaPsi,
            useLM = TRUE
        ),
        fullAE  = fitFraserAE(
            fds = fds,
            q = q,
            type = type,
            noiseAlpha = noiseAlpha,
            rhoRange = rhoRange,
            lambda = 0,
            convergence = convergence,
            iterations = iterations,
            initialize = initialize,
            nSubset = nSubset,
            weighted = weighted,
            control = control,
            BPPARAM = BPPARAM,
            verbose = verbose,
            subset = FALSE,
            nrDecoderBatches = 1
        ),
        PCA         = fitPCA(
            fds = fds,
            q = q,
            psiType = type,
            rhoRange = rhoRange,
            noiseAlpha = NULL,
            BPPARAM = BPPARAM,
            subset = FALSE
        ),
        'PCA-regression' = fitPCA(
            fds = fds,
            q = q,
            psiType = type,
            rhoRange = rhoRange,
            noiseAlpha = noiseAlpha,
            BPPARAM = BPPARAM,
            subset = TRUE,
            nSubset = nSubset,
            minDeltaPsi = minDeltaPsi
        ),
        BB          = fitBB(fds = fds, psiType = type)
    )

    return(fds)

}

needsHyperOpt <- function(method){
    switch(method,
        "PCA-BB-full"            = TRUE,
        "PCA-reg-full"           = TRUE,
        FraseR                   = TRUE,
        "PCA-BB-Decoder"         = TRUE,
        "PCA-BB-Decoder-no-weights" = TRUE,
        "FraseR-5DecoderBatches" = TRUE,
        "FraseR-1DecoderBatches" = TRUE,
        "FraseR-weighted"        = TRUE,
        fullFraseR               = TRUE,
        PCA                      = TRUE,
        'PCA-regression'         = TRUE,
        BB                       = FALSE,
        stop("Method not found: '", method, "'!")
    )
}



#'
#' Setting the hyper parameter optimization algorithm
#' for a given correction method
#'
#' @noRd
getHyperOptimCorrectionMethod <- function(correction){
    switch(correction,
            "PCA-BB-full"            = "PCA",
            "PCA-reg-full"           = "PCA",
            "PCA-BB-Decoder"         = "PCA-BB-Decoder",
            correction
    )
}

fitPCA <- function(fds, q, psiType, rhoRange=c(1e-5, 1-1e-5), noiseAlpha=NULL,
                    BPPARAM=bpparam(), subset=FALSE, minDeltaPsi,
                    nSubset=15000, useLM=FALSE){
    counts(fds, type=psiType, side="other", HDF5=FALSE) <- as.matrix(
            counts(fds, type=psiType, side="other"))
    counts(fds, type=psiType, side="ofInterest", HDF5=FALSE) <- as.matrix(
            counts(fds, type=psiType, side="ofInterest"))

    #+ subset fitting
    currentType(fds) <- psiType
    curDims <- dim(K(fds, psiType))
    if(!is.null(noiseAlpha)){
        noise(fds, type=type) <- matrix(rnorm(prod(curDims)), nrow=curDims[2])
    }

    # subset to fit the encoder
    fds_pca <- fds
    if(isTRUE(subset)){
        exMask <- variableJunctions(fds, type, minDeltaPsi=minDeltaPsi)
        fds_pca <- fds[exMask,,by=type]
        exMask2 <- subsetKMostVariableJunctions(fds_pca, type, nSubset)
        fds_pca <- fds_pca[exMask2,,by=type]
        featureExclusionMask(fds_pca) <- TRUE

        # set correct exclusion mask for x computation
        exMask[exMask == TRUE] <- exMask2
        featureExclusionMask(fds) <- exMask
    }

    # PCA on subset -> E matrix
    message(date(), ": Computing PCA ...")
    xin <- x(fds_pca, noiseAlpha=noiseAlpha, center=TRUE)
    pca <- pca(xin, nPcs=q)
    pc <- pcaMethods::loadings(pca)
    E(fds) <- pc

    # D and b on full matrix
    x <- x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)
    if(isTRUE(subset) | isTRUE(useLM)){
        # linear regression to fit D matrix
        lmFit  <- lm(x ~ H(fds))
        D(fds) <- t(lmFit$coefficients[-1,])
        b(fds) <- lmFit$coefficients[1,]
    } else{
        D(fds) <- pc
        b(fds) <- colMeans2(x)
    }

    # fit rho
    message(date(), ": Fitting rho ...")
    fds <- updateRho(fds, type=psiType, rhoRange=rhoRange,
            BPPARAM=BPPARAM, verbose=TRUE)

    metadata(fds)[[paste0('loss_', psiType)]] <- lossED(
            fds, byRows=TRUE, noiseAlpha=noiseAlpha)
    # store corrected logit psi
    predictedMeans(fds, psiType, withDimnames=FALSE) <- t(predictMu(fds))

    return(fds)
}

fitBB <- function(fds, psiType){
    currentType(fds) <- psiType
    fds <- pvalueByBetaBinomialPerType(fds=fds,
                                        aname=paste0("pvalues_BB_", psiType),
                                        psiType=psiType, 
                                        pvalFun=betabinVglmTest)
    # predictedMeans(fds, type=psiType) <- rowMeans(
    #         getAssayMatrix(fds, type=psiType))
    predictedMeans(fds, type=psiType) <-
        mcols(fds, type=psiType)[,paste0(psiType, "_alpha")] /
        ( mcols(fds, type=psiType)[,paste0(psiType, "_alpha")] +
                mcols(fds, type=psiType)[,paste0(psiType, "_beta")] )
    fds
}

fitFraserAE <- function(fds, q, type, noiseAlpha, rhoRange, lambda, convergence,
                        iterations, initialize, control, BPPARAM, verbose,
                        subset=TRUE, nrDecoderBatches, nSubset=15000, weighted,
                        latentSpace = 'AE'){
    #+ subset fitting
    curDims <- dim(K(fds, type))
    if(isTRUE(subset)){
        probE <- max(0.001, min(1,30000/curDims[1]))
        featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
                                            replace=TRUE, 
                                            prob=c(probE, 1-probE))
    } else{
        featureExclusionMask(fds) <- rep(TRUE,curDims[1])
    }
    print(table(featureExclusionMask(fds)))

    fds <- fitAutoencoder(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha,
                            rhoRange=rhoRange, lambda=lambda, 
                            convergence=convergence, iterations=iterations, 
                            initialize=initialize, nSubset=nSubset,
                            weighted=weighted, control=control, 
                            BPPARAM=BPPARAM, verbose=verbose, 
                            nrDecoderBatches=nrDecoderBatches, 
                            latentSpace=latentSpace )
    return(fds)

}
