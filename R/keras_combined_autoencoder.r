if(FALSE){
get_combined_keras_model <- function(inputDim, q, b, rho, optimizer, weights=NULL){
    if(FALSE){
        inputDim <- ncol(kt)
        inputDim
        q <- 10
        b <- bIn
        length(b)
        rho <- rho(fds)
        optimizer <- optimizer_adam(clipvalue=1, lr=lr)
        weights <- NULL
    }

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

    transformedData <- list(
            fds_inputK %>% layer_pseudo_k,
            fds_inputN %>% layer_pseudo_n) %>%
        layer_logit %>%
        layer_centered %>%
        list(inputNoise) %>%
        layer_noise

    #
    # mu model
    #
    # denoising autoencoder
    cmm <- rl2 <- NULL
    rl2     <- regularizer_l2()
    #cmm     <- constraint_minmaxnorm(-10, 10)
    encoder <- layer_dense(units=q,        kernel_regularizer=rl2, name="encoder", use_bias=FALSE)
    decoder <- layer_dense(units=inputDim, kernel_regularizer=rl2, name="decoder")

    # build model
    output_mu <- transformedData %>% encoder %>% decoder

    #
    # rho model
    #
    encoder_rho <- layer_dense(units=q,        kernel_regularizer=rl2, kernel_constraint=cmm, name="encoder_rho")
    decoder_rho <- layer_dense(units=inputDim, kernel_regularizer=rl2, kernel_constraint=cmm, name="decoder_rho")
    tf.math <- reticulate::import("tensorflow.math")
    layer_rho_collapse <- layer_lambda(name="rho_collapse", output_shape=c(1L, inputDim), f=function(x) k_mean(x, axis=1L, keepdims=TRUE) * 0)

    # build model
    # output_rho  <- transformedData %>% encoder_rho %>% decoder_rho %>% layer_rho_collapse %>%
    #    layer_lambda(name="rho_si", f=function(x) k_sigmoid(x))
    # output_rho <- k_variable(numeric(inputDim)) %>%
    output_rho  <- transformedData %>% decoder_rho %>% layer_rho_collapse %>%
        layer_dense(units=inputDim, name="rho_si", bias_constraint=constraint_minmaxnorm(1e-5, 1-1e-5, rate=1))
            # layer_lambda(name="rho_si", f=function(x) k_sigmoid(x)) + layer_add()

    # loss layer
    layer_loss <- layer_lambda(name="bb_rho_loss", f=function(x){
        keras_nll_betaBin(NULL, x[[1]], x[[2]], x[[3]], fds_rhoIn, mean=FALSE)
    })

    # output layer
    # output <- layer_concatenate(list(output_rho, output_mu), axis=0L)
    output <- list(output_mu, fds_inputK, fds_inputN, output_rho) %>% layer_loss
    output <- list(output_mu, fds_inputK, fds_inputN) %>% layer_loss

    # init bias for decoder
    set_weights(decoder, list(get_weights(decoder)[[1]], array(b)))
    # set_weights(decoder_rho, list(get_weights(decoder_rho)[[1]], array(rho)))

    # create and compile model
    fds_model <- keras_model(
        inputs = list(fds_inputK, fds_inputN, inputNoise),
        outputs = list(output))


    keras_loss_bb_combined <<- custom_metric("betabin_combined", function(y_true, y_pred){ k_mean(y_pred) })

    fds_model <- fds_model %>% compile(
        optimizer = optimizer,
        metrics = 'mean_squared_error',
        loss = keras_loss_bb_combined)
    model <- fds_model
    model
    # use known weights if provided
    if(!is.null(weights)){
        set_weigths(model, weights)
    }

    return(fds_model)
}

keras_nll_betaBin <- function(y_true, y_pred, k, n, rho, mean=TRUE){
    rho    <- k_clip(rho, 1e-5, 1-1e-5)
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

    if(isTRUE(mean)){
        return(k_mean(nllFull))
    }
    return(nllFull)
}


keras_loss_bb_combined <- custom_metric("betabin_combined", function(y_true, y_pred){
    keras_nll_betaBin(y_true[[1]], y_pred[[1]], k=fds_inputK, n=fds_inputN,
            rho=tf.math$reduce_mean(y_pred[[2]], axis=0L, keepdims=TRUE)) })

keras_loss_bb_combined <- custom_metric("betabin_combined", function(y_true, y_pred){
    keras_nll_betaBin(y_true[[0]], y_pred[[0]], k=fds_inputK, n=fds_inputN,
                      rho=y_pred[[1]]) })


keras_loss_bb_combined <- custom_metric("betabin_combined", function(y_true, y_pred){
    print(y_true[[1]])
    print(y_pred[[1]])
})

fds_model <- keras_model(
    inputs = list(fds_inputK, fds_inputN, inputNoise),
    outputs = list(output_mu, output_rho))
fds_model <- fds_model %>% compile(
    optimizer = optimizer,
    metrics = 'mean_squared_error',
    loss = keras_loss_bb_combined)
model <- fds_model


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
fit_keras_bb_dea_combined <- function(fds, q, type, noiseAlpha, rhoRange, BPPARAM, lr=0.0001, epochs=50, patience=5, reUseWeights=TRUE){
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
    rho(fds) <- methodOfMomentsRho(t(kt), t(nt))

    # get model
    message(date(), ": Init model ...")
    optimizer <- optimizer_adam(clipvalue=1, lr=lr)
    model <- get_keras_model(inputDim=ncol(kt), q=q, b=bIn, rho=rho(fds), optimizer=optimizer)

    # keep current weight setup
    curWeights  <- get_weights(model)
    initWeights <- get_weights(model)
    for(i in 1:10){

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
            callback_early_stopping(monitor="loss", patience=patience))
        history <- model %>% keras::fit(
            x = list(kt, nt, noise),
            #y = rbind(t(rho(fds)), xOut),
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

    return(fds)
}

get_prediction_model <- function(model){
    intermediate_layer_model = keras_model(
            input=list(
                get_layer(model, "k")$input,
                get_layer(model, "n")$input,
                get_layer(model, "noise")$input),
            output=list(
                get_layer(model, "decoder")$output)
                #,
                #get_layer(model, "rho_collapse")$output)
            )
    pred_mu <- intermediate_layer_model %>% predict(list(kt, nt, noise))
    length(pred_res)
    dim(pred_res[[1]])
    dim(pred_res[[2]])
    pred_res[[1]][1:10, 1:10]
    pred_res[[2]][1:10, 1:10]
}
}
