#'
#' ### MODEL
#'
fit_autoenc <- function(fds, type=currentType(fds), q_guess=round(ncol(fds)/4),
                    noiseRatio=0.5, BPPARAM=bpparam(), iterations=3){

    message(paste(date(), "; q:", q_guess, "; noise: ", noiseRatio))

    # setup object
    currentType(fds) <- type

    # train AE
    fds <- fitAutoencoder(fds, type=type, q=q_guess, iterations=iterations,
            verbose=TRUE, BPPARAM=BPPARAM)
    curLoss <- metadata(fds)[['loss']][length(metadata(fds)[['loss']])]

    return(list(fds=fds, evaluation=curLoss))
}

predict_outliers <- function(fds, type, BPPARAM){

    fds <- calculatePvalues(fds, type=type, BPPARAM=BPPARAM)

    return(fds)
}

eval_prot <- function(fds, type){
    index <- getSiteIndex(fds, type)
    idx   <- !duplicated(index)

    scores <- -as.vector(pVals(fds, type=type)[idx,])

    dt <- cbind(data.table(id=index), as.data.table(assay(fds, paste0("trueOutliers_", type))))
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


findEncodingDim <- function(i, fds, type, params, internalBPPARAM=SerialParam(), iterations){

    q_guess    <- params[i, "q"]
    noiseRatio <- params[i, "noise"]
    message(paste(i, ";\t", q_guess, ";\t", noiseRatio))

    res_fit <- fit_autoenc(fds=fds, type=type, q_guess=q_guess,
            noiseRatio=noiseRatio, BPPARAM=internalBPPARAM, iterations=iterations)
    res_pvals <- predict_outliers(res_fit$fds, type=type, BPPARAM=internalBPPARAM)
    evals <- eval_prot(res_pvals, type=type)

    return(list(q=q_guess, noiseRatio=noiseRatio, loss=res_fit$evaluation, aroc=evals))
}

optimHyperParams <- function(fds, type, q_param=seq(3, min(30, ncol(fds)), by=3),
                    noise_param=c(0, 0.5, 1, 2), minDeltaPsi=0.1, BPPARAM=bpparam(),
                    iterations=5, setSubset=0.1, injectFreq=1e-2){
    # make copy to return it later
    fds_copy <- fds
    currentType(fds) <- type

    #
    # put the most important stuff into memory
    #
    dontWriteHDF5(fds) <- TRUE
    counts(fds, type=type, side="other")      <- as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest") <- as.matrix(counts(fds, type=type, side="ofInterest"))

    #'
    #' remove non variable junctions
    #'
    psi <- K(fds)/N(fds)
    j2keep <- rowMaxs(abs(psi - rowMeans(psi, na.rm=TRUE)), na.rm=TRUE) >= minDeltaPsi
    message("dPsi filter:", paste(names(table(j2keep)), table(j2keep), collapse="\t", sep=": "))
    fds <- fds[j2keep,,by=type]

    #
    # subset for finding encoding dimensions
    #
    if(isScalarNumeric(setSubset) & nrow(fds) >= 500){
        exVec <- sample(c(T, F), nrow(mcols(fds, type=type)), replace=TRUE, prob=c(setSubset, 1-setSubset))
        featureExclusionMask(fds, type=type) <- exVec
        message("Exclusion matrix: ", paste(names(table(exVec)), table(exVec), collapse="\t", sep=": "))
    }
    fds <- fds[featureExclusionMask(fds, type=type),,by=type]

    # inject outliers
    fds <- injectOutliers(fds, type=type, freq=injectFreq, minDpsi=0.1, method="samplePSI")

    # run hyper parameter optimization
    params <- expand.grid(q=q_param, noise=noise_param)
    res <- bplapply(seq_len(nrow(params)), findEncodingDim, fds=fds, type=type,
                    iterations=iterations, params=params, BPPARAM=BPPARAM)

    # res2 <- res
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

