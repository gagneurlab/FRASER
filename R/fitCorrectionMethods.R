
fit <- function(fds, correction=c("FraseR", "fullFraseR", "PCA", "fullPCA", "PEER", "PEERdecoder", "BB", "kerasDAE", "kerasBBdAE"),
                    q, type="psi3", rhoRange=c(1e-5, 1-1e-5),
                    noiseAlpha=1, lambda=0, convergence=1e-5, iterations=15,
                    initialize=TRUE, control=list(), BPPARAM=bpparam(),
                    verbose=FALSE, recommendedPEERq=TRUE, lr=0.00005, epochs=20){
    method <- match.arg(correction)

    fds <- switch(method,
            FraseR      = fitFraserAE(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha, rhoRange=rhoRange, lambda=lambda,
                               convergence=convergence, iterations=iterations, initialize=initialize,
                               control=control, BPPARAM=BPPARAM, verbose=verbose, subset=TRUE),
            fullFraseR  = fitFraserAE(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha, rhoRange=rhoRange, lambda=lambda,
                                convergence=convergence, iterations=iterations, initialize=initialize,
                                control=control, BPPARAM=BPPARAM, verbose=verbose, subset=FALSE),
            PCA         = fitPCA(fds=fds, q=q, psiType=type, rhoRange=rhoRange, BPPARAM=BPPARAM, subset=TRUE),
            fullPCA     = fitPCA(fds=fds, q=q, psiType=type, rhoRange=rhoRange, BPPARAM=BPPARAM, subset=FALSE),
            PEER        = fitPEER(fds=fds, q=q, psiType=type, recomendedQ=recommendedPEERq, rhoRange=rhoRange, BPPARAM=BPPARAM),
            PEERdecoder = fitPEERDecoder(fds=fds, q=q, psiType=type, recomendedQ=recommendedPEERq, rhoRange=rhoRange, BPPARAM=BPPARAM),
            BB          = fitBB(fds=fds, psiType=type),
            kerasDAE    = fitKerasDAE(fds=fds, psiType=type, q=q, noiseAlpha=noiseAlpha, rhoRange=rhoRange, BPPARAM=BPPARAM),
            kerasBBdAE  = fit_keras_bb_dea(fds=fds, type=type, q=q, noiseAlpha=noiseAlpha,
                    rhoRange=rhoRange, BPPARAM=BPPARAM, lr=lr, patience=3, epochs=epochs, reUseWeights=FALSE, iterations=iterations))

  return(fds)

}

fitPCA <- function(fds, q, psiType, rhoRange=c(1e-5, 1-1e-5), BPPARAM=parallel(fds), subset=FALSE){

  #+ subset fitting
  currentType(fds) <- psiType
  curDims <- dim(K(fds, psiType))
  if(isTRUE(subset)){
    probE <- max(0.001, min(1,30000/curDims[1]))
    featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
                                        replace=TRUE, prob=c(probE, 1-probE))
  } else{
    featureExclusionMask(fds) <- rep(TRUE,curDims[1])
  }
  print(table(featureExclusionMask(fds)))

  # PCA on subset -> E matrix
  currentNoiseAlpha(fds) <- NULL
  message(date(), " Computing PCA ...")
  pca <- pca(as.matrix(x(fds, noiseAlpha=NULL, center=FALSE)), nPcs=q)
  pc <- pcaMethods::loadings(pca)
  E(fds) <- pc
  if(isTRUE(subset)){
    # linear regression to fit D matrix
    lmFit <- lm(as.matrix(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)) ~ as.matrix(H(fds)))
    D(fds) <- t(lmFit$coefficients[-1,])
    b(fds) <- lmFit$coefficients[1,]
  } else{
    D(fds) <- pc
    b(fds) <- double(nrow(mcols(fds, psiType)))
  }

  # fit rho
  message(date(), " Fitting rho ...")
  counts(fds, type=psiType, side="other", HDF5=FALSE) <- as.matrix(N(fds) - K(fds))
  counts(fds, type=psiType, side="ofInterest", HDF5=FALSE) <- as.matrix(K(fds))
  fds <- updateRho(fds, type=psiType, rhoRange=rhoRange, BPPARAM=BPPARAM, verbose=TRUE)
  # store corrected logit psi
  predictedMeans(fds, psiType) <- t(predictMu(fds))

  return(fds)
}

fitPEER <-function(fds, q, psiType, recomendedQ=TRUE, rhoRange=c(1e-5, 1-1e-5), BPPARAM=parallel(fds)){

  # set featureExclusionMask of all junctions to TRUE for peer
  currentType(fds) <- psiType
  featureExclusionMask(fds, type=psiType) <- rep(TRUE, nrow(mcols(fds, type=psiType)))

  #+ PEER
  require(peer)
  #+ prepare PEER model
  if(isTRUE(recomendedQ)){
    # recommendation by PEER: min(0.25*n, 100)
    q <- min(as.integer(0.25* ncol(fds)), 100)
  }
  maxFactors <- q                        # number of known hidden factors
  model <- PEER()
  PEER_setPhenoMean(model, as.matrix(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)))
  PEER_setNk(model, maxFactors)          # nr of hidden confounders
  PEER_setNmax_iterations(model, 1000)   # 1000 iterations is default
  # PEER_setAdd_mean(model, TRUE)        # should mean expression be added as additional factor? currently FALSE

  #+ run full Peer pipeline
  message(date(), "Fitting PEER model ...")
  PEER_update(model)

  #+ extract PEER data
  peerResiduals <- PEER_getResiduals(model)
  peerLogitMu <- t(as.matrix(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)) - peerResiduals)

  #+ save peer model in fds object
  setAssayMatrix(fds, "peerLogitMu", type=psiType) <- peerLogitMu
  predictedMeans(fds, psiType) <- predictMuCpp(peerLogitMu)
  metadata(fds)[[paste0("PEERmodel_", psiType)]] <- list(
    alpha     = PEER_getAlpha(model),
    residuals = PEER_getResiduals(model),
    W         = PEER_getW(model),
    hiddenSpace = PEER_getX(model))

  #+ fit rho
  message(date(), "Fitting rho ...")
  k <- as.matrix(K(fds, psiType))
  n <- as.matrix(N(fds, psiType))
  y <- peerLogitMu
  fitparameters <- bplapply(seq_len(nrow(k)), estRho,
                            k=k, n=n, y=y, rhoRange=rhoRange,
                            BPPARAM=BPPARAM, nll=truncNLL_rho)
  rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
  print(summary(rho(fds)))

  return(fds)

}

fitPEERDecoder <-function(fds, q, psiType, recomendedQ=TRUE, rhoRange=c(1e-5, 1-1e-5), BPPARAM=parallel(fds)){

  # set featureExclusionMask of all junctions to TRUE for peer
  currentType(fds) <- psiType
  featureExclusionMask(fds, type=psiType) <- rep(TRUE, nrow(mcols(fds, type=psiType)))

  #+ PEER
  require(peer)
  #+ prepare PEER model
  maxFactors <- q                        # number of known hidden factors
  if(isTRUE(recomendedQ)){
    # recommendation by PEER: min(0.25*n, 100)
    maxFactors <- min(as.integer(0.25* ncol(fds)), 100)
  }
  model <- PEER()
  x <- x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)
  PEER_setPhenoMean(model, as.matrix(x))
  PEER_setNk(model, maxFactors)          # nr of hidden confounders
  PEER_setNmax_iterations(model, 1000)   # 1000 iterations is default
  # PEER_setAdd_mean(model, TRUE)        # should mean expression be added as additional factor? currently FALSE

  #+ run full Peer pipeline
  message(date(), "Fitting PEER model ...")
  PEER_update(model)

  #+ extract PEER data
  peerResiduals <- PEER_getResiduals(model)
  peerLogitMu <- t(as.matrix(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)) - peerResiduals)

  #+ save peer model in fds object
  setAssayMatrix(fds, "peerLogitMu", type=psiType) <- peerLogitMu
  metadata(fds)[[paste0("PEERmodel_", psiType)]] <- list(
    alpha     = PEER_getAlpha(model),
    residuals = PEER_getResiduals(model),
    W         = PEER_getW(model),
    hiddenSpace = PEER_getX(model))

  #+ fit rho
  message(date(), "Fitting rho ...")
  k <- as.matrix(K(fds, psiType))
  n <- as.matrix(N(fds, psiType))
  y <- peerLogitMu
  fitparameters <- bplapply(seq_len(nrow(k)), estRho,
                            k=k, n=n, y=y, rhoRange=rhoRange,
                            BPPARAM=BPPARAM, nll=truncNLL_rho)
  rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
  print(summary(rho(fds)))

  #+ fit decoder on peer hidden space
  D <- matrix(rnorm(prod(dim(fds))), nrow=nrow(fds), ncol=q)
  b <- colMeans2(x)
  rho <- rho(fds)
  alphas <- PEER_getAlpha(model)[,1]
  H <- PEER_getX(model)[,order(alphas, decreasing = TRUE)[1:q]]

  # run fits
  fitls <- bplapply(seq_len(nrow(k)), singleDFit, D=D, b=b, k=k, n=n, H=H,
                    rho=rho, lambda=0, control=list(), BPPARAM=parallel(fds))

  # extract infos
  parMat <- vapply(fitls, '[[', double(ncol(D) + 1), 'par')

  # update b and D
  b(fds) <- parMat[1,]
  D(fds) <- t(parMat)[,-1]

  # predict and save (logit) psi
  y <- predictYCpp(as.matrix(H), D(fds), b(fds))
  predictedMeans(fds, psiType) <- t(predictMuCpp(y))

  return(fds)

}

fitBB <- function(fds, psiType){
  currentType(fds) <- psiType
  fds <- pvalueByBetaBinomialPerType(fds=fds, aname=paste0("pvalue_", psiType), psiType=psiType, pvalFun=betabinVglmTest)

  return(fds)
}

fitFraserAE <- function(fds, q, type, noiseAlpha, rhoRange, lambda, convergence, iterations, initialize,
                                  control, BPPARAM, verbose, subset=TRUE){
  #+ subset fitting
  curDims <- dim(K(fds, type))
  if(isTRUE(subset)){
    probE <- max(0.001, min(1,30000/curDims[1]))
    featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
                                        replace=TRUE, prob=c(probE, 1-probE))
  } else{
    featureExclusionMask(fds) <- rep(TRUE,curDims[1])
  }
  print(table(featureExclusionMask(fds)))

  fds <- fitAutoencoder(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha, rhoRange=rhoRange, lambda=lambda,
                 convergence=convergence, iterations=iterations, initialize=initialize,
                 control=control, BPPARAM=BPPARAM, verbose=verbose)
  return(fds)

}


fitKerasDAE <-function(fds, q, psiType, rhoRange=c(1e-5, 1-1e-5), noiseAlpha=1, BPPARAM=parallel(fds),
                    pypath="/opt/modules/i12g/anaconda/3-5.0.1/envs/omicsOUTRIDER/bin/python"){

    # set needed default values
    message(date(), "set params for kerasDAE ...")
    currentType(fds) <- psiType
    featureExclusionMask(fds, type=psiType) <- rep(TRUE, nrow(mcols(fds, type=psiType)))
    currentNoiseAlpha(fds) <- noiseAlpha
    noise(fds) <- matrix(rnorm(prod(dim(t(K(fds))))), ncol=ncol(t(K(fds))))

    # get data into memory
    counts(fds, type=type, side="other", HDF5=FALSE)      <- as.matrix(N(fds) - K(fds))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(K(fds))

    #+ setup keras and python
    require(keras)
    use_python(pypath)

    # data keras input
    X      <- x(fds, noiseAlpha=NULL, all=TRUE, center=FALSE)
    X_corr <- x(fds, noiseAlpha=currentNoiseAlpha(fds), all=TRUE)
    bIn    <- colMeans(X)

    # shuffel input
    sidx <- sample(nrow(X), nrow(X))
    X <- X[sidx,]
    X_corr <- X_corr[sidx,]

    # layers
    message(date(), "setup keras model ...")
    input   <- layer_input(ncol(X))
    encoder <- layer_dense(units=q)
    decoder <- layer_dense(units=ncol(X))

    # full model
    output <- input %>% encoder %>% decoder

    # init bias for decoder
    set_weights(decoder, list(get_weights(decoder)[[1]], array(bIn)))

    # create and compile model
    model <- keras_model(
            inputs = list(input),
            outputs = output) %>%
        compile(
            optimizer = 'rmsprop',
            loss = 'mean_squared_error',
            metrics = c('mean_squared_error'))

    # fit the model on whole data with early stopping
    message(date(), "fitting keras model ...")
    cb_es <- list(callback_early_stopping(monitor="val_loss", patience=10))
    history <- model %>% keras::fit(x=X_corr, y=X,
        epochs=300, batch_size=16,
        callbacks=cb_es, validation_split=0.2)
    print(plot(history, main = "History of fitting the model, batch-aware PCA", method="ggplot2") +
        scale_y_log10() +
        xlim(0, cb_es[[1]]$stopped_epoch))

    # get predictions
    pred_mu <- model %>% predict(X_corr)

    # plot predictions vs observerd
    d2p <- data.table(x=as.vector(X), mu=as.vector(pred_mu))[sample(.N, 1e5)]
    require(LSD)
    heatscatter(d2p[,mu], d2p[,x], ylab="observed", xlab="predicted",
            main="predictions vs observed\nrandom 100k points")
    grid()
    abline(0,1,col="red")

    # set predictions
    predictedMeans(fds, psiType) <- predictMuCpp(t(pred_mu))

    # fit rho
    message(date(), "fitting rho ...")
    fitparameters <- bplapply(seq_len(nrow(K(fds))),
            estRho, k=K(fds), n=N(fds), y=t(pred_mu),
            rhoRange=rhoRange, BPPARAM=BPPARAM,
            nll=truncNLL_rho)
    rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
    print(summary(rho(fds)))

    return(fds)
}
