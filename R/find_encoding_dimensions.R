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
                    internalBPPARAM=1, iterations,
                    nrDecoderBatches){
    iBPPARAM <- MulticoreParam(internalBPPARAM)

    q_guess    <- params[i, "q"]
    noiseRatio <- params[i, "noise"]
    message(paste(i, ";\t", q_guess, ";\t", noiseRatio))
    correction <- getHyperOptimCorrectionMethod(correction)

    res_fit <- fit_autoenc(fds=fds, type=type, q_guess=q_guess,
            correction=correction, nrDecoderBatches=nrDecoderBatches,
            noiseRatio=noiseRatio, BPPARAM=iBPPARAM,
            iterations=iterations)
    res_pvals <- predict_outliers(res_fit$fds, correction=correction,
            type=type, BPPARAM=iBPPARAM)
    evals <- eval_prot(res_pvals, type=type)

    return(list(q=q_guess, noiseRatio=noiseRatio, loss=res_fit$evaluation, aroc=evals))
}

optimHyperParams <- function(fds, type, correction,
                    q_param=seq(2, min(40, ncol(fds)), by=3),
                    noise_param=c(0, 0.5, 1, 2, 5), minDeltaPsi=0.1,
                    iterations=5, setSubset=15000, injectFreq=1e-2,
                    nrDecoderBatches=1, BPPARAM=bpparam(), internalThreads=3){
    if(isFALSE(needsHyperOpt(correction))){
        message(date(), ": For correction '", correction, "' no hyper paramter",
                "optimization is needed.")
        data <- data.table(q=NA, noise=0, eval=1, aroc=1)
        hyperParams(fds, type=type) <- data
        return(fds)
    }

    #
    # put the most important stuff into memory
    #
    currentType(fds) <- type
    counts(fds, type=type, side="other", HDF5=FALSE) <-
            as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <-
            as.matrix(counts(fds, type=type, side="ofInterest"))

    #'
    #' remove non variable and low abundance junctions
    #'
    j2keepVa <- variableJunctions(fds, type, minDeltaPsi)
    j2keepDP <- rowQuantiles(K(fds, type), probs=0.75) >= 10
    j2keep <- j2keepDP & j2keepVa
    message("dPsi filter:", pasteTable(j2keep))
    # TODO fds <- fds[j2keep,,by=type]

    optData <- data.table()
    for(nsub in setSubset){
        #
        # subset for finding encoding dimensions
        # most variable functions + random subset for decoder
        exMask <- subsetKMostVariableJunctions(fds[j2keep,,by=type], type, nsub)
        j2keep[j2keep==TRUE] <- exMask

        # keep n most variable junctions + random subset
        j2keep <- j2keep | sample(c(TRUE, FALSE), length(j2keep), replace=TRUE,
                prob=c(nsub/length(j2keep), 1 - nsub/length(j2keep)))
        message("Exclusion matrix: ", pasteTable(j2keep))

        # make copy for testing
        fds_copy <- fds
        dontWriteHDF5(fds_copy) <- TRUE
        featureExclusionMask(fds_copy) <- j2keep
        fds_copy <- fds_copy[j2keep,,by=type]
        currentType(fds_copy) <- type

        # inject outliers
        fds_copy <- injectOutliers(fds_copy, type=type, freq=injectFreq,
                minDpsi=minDeltaPsi, method="samplePSI", BPPARAM=BPPARAM)

        # remove unneeded blocks to save memory
        a2rm <- paste(sep="_", c("originalCounts", "originalOtherCounts"),
                rep(psiTypes, 2))
        for(a in a2rm){
            assay(fds_copy, a) <- NULL
        }
        metadata(fds_copy) <- list()
        gc()

        # reset lost important values
        currentType(fds_copy) <- type
        dontWriteHDF5(fds_copy) <- TRUE

        # run hyper parameter optimization
        params <- expand.grid(q=q_param, noise=noise_param)
        message(date(), ": Run hyper optimization with ", nrow(params), " options.")
        res <- bplapply(seq_len(nrow(params)), findEncodingDim, fds=fds_copy, type=type,
                        iterations=iterations, params=params, correction=correction,
                        BPPARAM=BPPARAM, nrDecoderBatches=nrDecoderBatches,
                        internalBPPARAM=internalThreads)

        data <- data.table(
            q=sapply(res, "[[", "q"),
            noise=sapply(res, "[[", "noiseRatio"),
            nsubset=nsub,
            eval=sapply(res, "[[", "loss"),
            aroc=sapply(res, "[[", "aroc"))

        optData <- rbind(optData, data)
    }

    print(plot_find_enc_results(optData))

    hyperParams(fds, type=type) <- optData
    return(fds)
}

plot_find_enc_results <- function(data){
    if(!"nsubset" %in% colnames(data)){
        data[,nsubset:="NA"]
    }
    data[,noise:=as.factor(noise)]
    data[,nsubset:=as.factor(nsubset)]

    g1 <- ggplot(data, aes(q, aroc, col=nsubset, linetype=noise)) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        ggtitle("Q estimation") +
        xlab("Estimated q") +
        ylab("Area under the ROC curve")

    g2 <- ggplot(data, aes(q, eval, col=nsubset, linetype=noise)) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        ggtitle("Q estimation") +
        xlab("Estimated q") +
        ylab("Model loss")


    ggarrange(g1, g2, nrow=2)
}

