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

findBestEncoding <- function(i, fds, type, params, iterations){

    BPPARAM <- SerialParam()
    q_guess    <- params[i, "q"]
    noiseRatio <- params[i, "noise"]
    message(paste(i, ";\t", q_guess, ";\t", noiseRatio))

    res_fit <- fit_autoenc(fds=fds, type=type, q_guess=q_guess,
            noiseRatio=noiseRatio, BPPARAM=BPPARAM, iterations=iterations)
    res_pvals <- predict_outliers(res_fit$fds, type=type, BPPARAM=BPPARAM)
    evals <- eval_prot(res_pvals, type=type)

    return(list(q=q_guess, noiseRatio=noiseRatio, loss=res_fit$evaluation, aroc=evals))
}

plot_find_enc_results <- function(data){
    require(ggplot2)
    require(patchwork)
    g1 <- ggplot(data, aes(q, aroc, col=as.factor(noise))) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        title("Q estimation") +
        xlab("Estimated q") +
        ylab("Area under the ROC curve")

    g2 <- ggplot(data, aes(q, eval, col=as.factor(noise))) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        title("Q estimation") +
        xlab("Estimated q") +
        ylab("Model loss")


    g1 + g2 + plot_layout(ncol=1, nrow=2)
}

