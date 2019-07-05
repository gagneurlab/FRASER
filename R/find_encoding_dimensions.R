#'
#' ### MODEL
#'
fit_autoenc <- function(fds, type=currentType(fds), q_guess=round(ncol(fds)/4),
                    correction, noiseRatio=0.5, BPPARAM=bpparam(),
                    iterations=3, nrDecoderBatches=1, verbose=FALSE){

    message(paste(date(), "; q:", q_guess, "; noise: ", noiseRatio))

    # setup object
    currentType(fds) <- type

    # train AE
    fds <- fit(fds, type=type, q=q_guess, correction=correction,
            iterations=iterations, nrDecoderBatches=nrDecoderBatches,
            verbose=verbose, BPPARAM=BPPARAM)
    curLoss <- metadata(fds)[[paste0('loss_', type)]]
    curLoss <- mean(curLoss[,ncol(curLoss)])

    return(list(fds=fds, evaluation=curLoss))
}

predict_outliers <- function(fds, type, correction, BPPARAM){

    fds <- calculatePvalues(fds, type=type, correction=correction,
            BPPARAM=BPPARAM)

    return(fds)
}

eval_prot <- function(fds, type){
    index <- getSiteIndex(fds, type)
    idx   <- !duplicated(index)

    scores <- -as.vector(pVals(fds, type=type)[idx,])

    dt <- cbind(data.table(id=index),
            as.data.table(assay(fds, paste0("trueOutliers_", type))))
    setkey(dt, id)
    labels <- as.vector(sapply(samples(fds), function(i){
        dttmp <- dt[,any(get(i) != 0),by=id]
        setkey(dttmp, id)
        dttmp[J(unique(index)), V1]
    })) + 0

    if(any(is.na(scores))){
        warning(sum(is.na(scores)), " P-values where NAs.")
        scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
    }
    pr <- pr.curve(scores, weights.class0=labels)
    pr

    return(pr$auc.integral)
}


findEncodingDim <- function(i, fds, type, params, correction,
                    internalBPPARAM=SerialParam(), iterations,
                    nrDecoderBatches){

    q_guess    <- params[i, "q"]
    noiseRatio <- params[i, "noise"]
    message(paste(i, ";\t", q_guess, ";\t", noiseRatio))

    res_fit <- fit_autoenc(fds=fds, type=type, q_guess=q_guess,
            correction=correction, nrDecoderBatches=nrDecoderBatches,
            noiseRatio=noiseRatio, BPPARAM=internalBPPARAM,
            iterations=iterations)
    res_pvals <- predict_outliers(res_fit$fds, correction=correction,
            type=type, BPPARAM=internalBPPARAM)
    evals <- eval_prot(res_pvals, type=type)

    return(list(q=q_guess, noiseRatio=noiseRatio, loss=res_fit$evaluation, aroc=evals))
}

optimHyperParams <- function(fds, type, correction,
                    q_param=seq(2, min(40, ncol(fds)), by=3),
                    noise_param=c(0, 0.5, 1, 2, 5), minDeltaPsi=0.1,
                    iterations=5, setSubset=15000, injectFreq=1e-2,
                    nrDecoderBatches=1, BPPARAM=bpparam()){
    if(needsHyperOpt(correction) == FALSE){
        message(date(), ": For correction '", correction, "' no hyper paramter",
                "optimization is needed.")
        data <- data.table(q=NA, noise=0, eval=1, aroc=1)
        hyperParams(fds, type=type) <- data
        return(fds)
    }

    # make copy to return it later
    fds_copy <- fds
    currentType(fds) <- type

    #
    # put the most important stuff into memory
    #
    dontWriteHDF5(fds) <- TRUE
    counts(fds, type=type, side="other")      <-
            as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest") <-
            as.matrix(counts(fds, type=type, side="ofInterest"))

    #'
    #' remove non variable and low abundance junctions
    #'
    j2keepVa <- variableJunctions(fds, type, minDeltaPsi)
    j2keepDP <- rowMedians(K(fds, type)) >= 5
    j2keep <- j2keepDP & j2keepVa
    message("dPsi filter:", pasteTable(j2keep))
    fds <- fds[j2keep,,by=type]

    #
    # subset for finding encoding dimensions
    #
    exMask <- subsetKMostVariableJunctions(fds, type, setSubset)
    featureExclusionMask(fds) <- exMask
    fds <- fds[exMask,,by=type]
    message("Exclusion matrix: ", pasteTable(exMask))

    # inject outliers
    fds <- injectOutliers(fds, type=type, freq=injectFreq,
            minDpsi=minDeltaPsi, method="samplePSI")

    # run hyper parameter optimization
    params <- expand.grid(q=q_param, noise=noise_param)
    res <- bplapply(seq_len(nrow(params)), findEncodingDim, fds=fds, type=type,
                    iterations=iterations, params=params, correction=correction,
                    BPPARAM=BPPARAM, nrDecoderBatches=nrDecoderBatches)

    data <- data.table(
        q=sapply(res, "[[", "q"),
        noise=sapply(res, "[[", "noiseRatio"),
        eval=sapply(res, "[[", "loss"),
        aroc=sapply(res, "[[", "aroc"))
    print(plot_find_enc_results(data))

    hyperParams(fds_copy, type=type) <- data
    return(fds_copy)
}

plot_find_enc_results <- function(data){
    g1 <- ggplot(data, aes(q, aroc, col=as.factor(noise))) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        ggtitle("Q estimation") +
        xlab("Estimated q") +
        ylab("Area under the ROC curve")

    g2 <- ggplot(data, aes(q, eval, col=as.factor(noise))) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        ggtitle("Q estimation") +
        xlab("Estimated q") +
        ylab("Model loss")


    g1 + g2 + plot_layout(ncol=1, nrow=2)
}

