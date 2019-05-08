get_keras_model <- function(inputDim, q, b, rho, optimizer, weights=NULL){

    # layers
    fds_inputK  <<- layer_input(inputDim, name="k")
    fds_inputN  <<- layer_input(inputDim, name="n")
    fds_rhoIn   <<- k_variable(rho,       name="rhoin")
    inputNoise   <- layer_input(inputDim, name="noise")
    bias         <- k_variable(b,         name="bias")

    # transformation layers
    layer_pseudo_n <- layer_lambda(name="pseudo_n",    f=function(x) x+2)
    layer_pseudo_k <- layer_lambda(name="pseudo_k",    f=function(x) x+1)
    layer_logit    <- layer_lambda(name="logit",       f=function(x) k_log((x[[1]]/x[[2]])/(1-(x[[1]]/x[[2]]))))
    layer_centered <- layer_lambda(name="rowcentered", f=function(x) x - bias)
    layer_noise    <- layer_lambda(name="corrupted",   f=function(x) x[[1]] + x[[2]])

    # denoising autoencoder
    rl2     <- regularizer_l2()
    cmm     <- constraint_minmaxnorm(-10, 10)
    encoder <- layer_dense(units=q,        kernel_regularizer=rl2, kernel_constraint=cmm, name="encoder")
    decoder <- layer_dense(units=inputDim, kernel_regularizer=rl2, kernel_constraint=cmm, name="decoder")

    # full model
    transformedData <- list(
        fds_inputK %>% layer_pseudo_k,
        fds_inputN %>% layer_pseudo_n) %>%
        layer_logit %>%
        layer_centered
    output <- list(transformedData, inputNoise) %>%
        layer_noise %>% encoder %>% decoder

    # init bias for decoder
    set_weights(decoder, list(get_weights(decoder)[[1]], array(b)))

    # create and compile model
    fds_model <- keras_model(
        inputs = list(fds_inputK, fds_inputN, inputNoise),
        outputs = output)
    fds_model <- fds_model %>% compile(
        optimizer = optimizer,
        metrics = 'mean_squared_error',
        loss = keras_loss_bb)

    # use known weights if provided
    if(!is.null(weights)){
        set_weigths(model, weights)
    }

    return(fds_model)
}

keras_nll_betaBin <- function(y_true, y_pred, k, n, rho){
    mu     <- k_exp(y_pred)/(1 + k_exp(y_pred))
    a      <- mu * (1/rho -1)
    b      <- ((mu - 1)*(rho - 1))/rho

    tf.lg <- tensorflow:::lgamma.tensorflow.tensor
    nllFull <- (
        - tf.lg(n + 3)
        + tf.lg(k + 2)
        + tf.lg(n - k + 2)

        - tf.lg(k + a + 1)
        - tf.lg(n - k + b + 1)
        + tf.lg(n + a + b + 2)

        - tf.lg(a + b)
        + tf.lg(a)
        + tf.lg(b)
    )

    k_mean(nllFull)
}

keras_loss_bb <- custom_metric("betabin", function(y_true, y_pred){
    keras_nll_betaBin(y_true, y_pred, k=fds_inputK, n=fds_inputN,
            rho=k_variable(fds_rhoIn)) })

#'
#'
#' Keras denoising autoencoder with betabinomial loss
#'
#' Install environment: (bash commands)
#'
#'   envName=omicsOUTRIDER
#'   module load i12g/anaconda/3-5.0.1
#'   conda create -n $envName python=3.6
#'   /opt/modules/i12g/anaconda/3-5.0.1/envs/$envName/bin/pip install -U keras tensorflow jupyterlab
#'
#' To use it run before training:
#'
#'   library(keras)
#'   use_python("/opt/modules/i12g/anaconda/3-5.0.1/envs/omicsOUTRIDER/bin/python", required=TRUE)
#'
fit_keras_bb_dea <- function(fds, q, type, noiseAlpha, rhoRange, BPPARAM, lr=0.0001, epochs=50, patience=5, reUseWeights=TRUE, iterations=10){
    currentType(fds)   <- type

    counts(fds, type=type, side="other", HDF5=FALSE)      <- as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(counts(fds, type=type, side="ofInterest"))

    # data keras input
    kt       <- t(K(fds))
    nt       <- t(N(fds))
    xMat     <- x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)
    xOut     <- x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)
    bIn      <- colMeans(xMat)
    noise    <- t(colSds(xMat) * noiseAlpha * t(matrix(rnorm(prod(dim(kt))), ncol=ncol(kt))))
    rho(fds) <- methodOfMomemtsRho(t(kt))

    # get model
    message(date(), ": Init model ...")
    optimizer <- optimizer_adam(clipvalue=1, lr=lr)
    model <- get_keras_model(inputDim=ncol(kt), q=q, b=bIn, rho=rho(fds), optimizer=optimizer)

    # keep current weight setup
    curWeights  <- get_weights(model)
    initWeights <- get_weights(model)
    for(i in seq_len(10)){

        # shuffel input
        sidx <- sample(nrow(kt), nrow(kt))
        kt <- kt[sidx,]
        nt <- nt[sidx,]

        # fit rho
        message(date(), ": Fit rho ...")
        pred_mu <- model %>% predict(list(kt, nt, noise))
        fitparameters <- bplapply(seq_len(ncol(kt)), estRho, nll=truncNLL_rho,
                k=t(kt), n=t(nt), y=t(pred_mu), rhoRange=rhoRange,
                BPPARAM=BPPARAM)
        rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
        k_set_value(fds_rhoIn, rho(fds))

        # train dAE
        message(date(), ": fit model ...")
        model <- get_keras_model(inputDim=ncol(kt), q=q, b=bIn,
                rho=rho(fds), optimizer=optimizer)
        if(isTRUE(reUseWeights)){
            set_weights(model) <- curWeights
        }
        cb_list <- list(
            callback_terminate_on_naan(),
            callback_early_stopping(patience=patience))
        history <- model %>% keras::fit(
            x = list(kt, nt, noise),
            y = xOut,
            shuffle = TRUE,
            epochs=epochs,
            batch_size=ceiling(nrow(kt)/10),
            validation_split=0.2,
            callbacks=cb_list)

        # save weights
        curWeights <- get_weights(model)

        message("Current loss:", rev(history$metrics$val_loss)[1])
        curEval <- model %>% evaluate(x=list(kt, nt, noise), y=xOut)
        message(paste(names(curEval), curEval, sep=":\t", collapse="\n"))
    }

    # predict mu and save it
    pred_mu <- model %>% predict(list(kt, nt, noise))
    predictedMeans(fds, type) <- predictMuCpp(t(pred_mu))

    # save weights
    metadata(fds)[[paste("dAE_keras_bb_weights_", type)]] <- curWeights
    
    # save noise
    noise(fds, type) <- noise

    return(fds)
}




