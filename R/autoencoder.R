#'
#' Main autoencoder fit function
#'
#' @noRd
#'
#' @export
fitAutoencoder <- function(fds, q, type="psi3", noiseAlpha=1, rhoRange=c(1e-5, 1-1e-5), lambda=0,
                           convergence=1e-5, iterations=15, initialize=TRUE,
                           control=list(), BPPARAM=bpparam(), verbose=FALSE){

    if(!bpisup(BPPARAM)){
        bpstart(BPPARAM)
    }
    dims <- c(row=nrow(mcols(fds, type=type)), col=ncol(fds))

    # set alpha for noise injection for denoising AE
    currentNoiseAlpha(fds) <- noiseAlpha
    curNoise <- rnorm(prod(dims), mean=0, sd=1) * noiseAlpha
    curNoise <- matrix(curNoise, nrow=dims[1], ncol=dims[2])
    noise(fds, type=type) <- curNoise

    # make sure its only in-memory data for k and n
    currentType(fds) <- type
    counts(fds, type=type, side="other", HDF5=FALSE) <- as.matrix(N(fds) - K(fds))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(K(fds))

    # copy fds object to save original input object
    # and create second object with only the subset to fit the encoder
    copy_fds <- fds
    fds <- fds[featureExclusionMask(fds, type=type),by=checkReadType(fds, type=type)]

    # subset noise so that it works with subsetted x
    noise(fds, type=type) <- curNoise[featureExclusionMask(copy_fds, type=type),]

    # initialize E and D using PCA and bias as zeros.
    if(isTRUE(initialize) | is.null(E(fds)) | is.null(D(fds))){
        fds <- initAutoencoder(fds, q, rhoRange, type=type)
    }

    # initial loss
    lossList <- c(init_pca=lossED(fds, lambda))
    print(paste0('Initial PCA loss: ', lossList[1]))

    #inital rho values
    if(isTRUE(verbose)){
        print(summary(rho(fds)))
    }

    # initialize D
    fds <- updateD(fds, type=type, lambda=lambda, control=control, BPPARAM=BPPARAM, verbose=verbose)
    lossList <- updateLossList(fds, lossList, 'init', 'D', lambda, verbose=verbose)

    # initialize rho step
    fds <- updateRho(fds, type=type, rhoRange, BPPARAM=BPPARAM, verbose=verbose)
    lossList <- updateLossList(fds, lossList, 'init', 'Rho', lambda, verbose=verbose)

    # optimize log likelihood
    t1 <- Sys.time()
    currentLoss <- lossED(fds, lambda)
    for(i in seq_len(iterations)){
        t2 <- Sys.time()

        # update E step
        fds <- updateE(fds, control=control, BPPARAM=BPPARAM, verbose=verbose)
        lossList <- updateLossList(fds, lossList, i, 'E', lambda, verbose=verbose)

        # update D step
        fds <- updateD(fds, type=type, lambda=lambda, control=control, BPPARAM=BPPARAM, verbose=verbose)
        lossList <- updateLossList(fds, lossList, i, 'D', lambda, verbose=verbose)

        # update rho step
        fds <- updateRho(fds, type=type, rhoRange, BPPARAM=BPPARAM, verbose=verbose)
        lossList <- updateLossList(fds, lossList, i, 'Rho', lambda, verbose=verbose)

        if(isTRUE(verbose)){
          print(paste('Time for one autoencoder loop:', Sys.time() - t2))
        } else {
          print(paste0(date(), ': Iteration: ', i, ' loss: ',
                       lossList[length(lossList)]))
        }

        # check
        curLossDiff <- abs(currentLoss - lossList[length(lossList) - 2:0])
        if(all(curLossDiff < convergence)){
          message(date(), ': the AE correction converged with:',
                  lossList[length(lossList)])
          break
        }
        currentLoss <- lossList[length(lossList)]
    }

    print(Sys.time() - t1)

    if(nrow(fds) == nrow(copy_fds)){    # TODO when using all features (!=nrow(fds)) in fitting for SE: also stop here and set copy_fds to fitted fds but with all junctions
                                        # (at the moment: fds contains only a subset of all junctions, so copy_fds <- fds doesn't work for SE)
        copy_fds <- fds
    } else {
        # update the D matrix and theta
        print("Finished with fitting the E matrix. Starting now with the D fit. ...")

        # set noiseAlpha to 0 or NULL to NOT use noise now (latent space already fitted)
        currentNoiseAlpha(fds) <- NULL

        copy_fds <- initAutoencoder(copy_fds, q, rhoRange, type=type)
        E(copy_fds) <- E(fds)

        for(i in seq_len(iterations)){
            t2 <- Sys.time()

            # update D step
            copy_fds <- updateD(copy_fds, type=type, lambda=lambda, control=control, BPPARAM=BPPARAM, verbose=verbose)
            lossList <- updateLossList(copy_fds, lossList, paste0("final_", i), 'D', lambda, verbose=verbose)

            # update rho step
            copy_fds <- updateRho(copy_fds, type=type, rhoRange, BPPARAM=BPPARAM, verbose=verbose)
            lossList <- updateLossList(copy_fds, lossList, paste0("final_", i), 'Rho', lambda, verbose=verbose)

            if(isTRUE(verbose)){
                print(paste('Time for one D & Rho loop:', Sys.time() - t2))
            } else {
                print(paste0(date(), ': Iteration: final_', i, ' loss: ',
                             lossList[length(lossList)]))
            }

            # check
            curLossDiff <- abs(currentLoss - lossList[length(lossList) - 1:0])
            if(all(curLossDiff < convergence)){
                message(date(), ': the final AE correction converged with:',
                        lossList[length(lossList)])
                break
            }
            currentLoss <- lossList[length(lossList)]
        }
    }

    print(paste0(i, ' Final betabin-AE loss: ', lossList[length(lossList)]))
    bpstop(BPPARAM)

    # add additional values for the user to the object
    metadata(copy_fds)[['dim']] <- dim(copy_fds)
    metadata(copy_fds)[['loss']] <- lossList
    metadata(copy_fds)[['convList']] <- lossList


    # add correction factors
    predictedMeans <- t(predictMu(copy_fds))
    stopifnot(identical(dim(K(copy_fds)), dim(predictedMeans)))
    predictedMeans(copy_fds, type) <- predictedMeans

    # validate object
    validObject(copy_fds)
    return(copy_fds)
}

initAutoencoder <- function(fds, q, rhoRange, type){

    pca <- pca(as.matrix(x(fds, all=TRUE)), nPcs=q)
    pc  <- pcaMethods::loadings(pca)
    #pc = matrix(rnorm(q*nrow(mcols(fds, type=type)), sd=0.1), ncol=q)

    # Set initial values from PCA
    D(fds) <- pc
    E(fds) <- pc[featureExclusionMask(fds),]
    b(fds) <- double(nrow(mcols(fds, type=type)))

    # initialize rho
    rho(fds) <- methodOfMomemtsRho(K(fds))

    # reset counters
    mcols(fds, type=type)[paste0('NumConvergedD_', type)] <- 0

    return(fds)
}

updateLossList <- function(fds, lossList, i, stepText, lambda, verbose){
    currLoss <- lossED(fds, lambda)
    lossList <- c(lossList, currLoss)
    names(lossList)[length(lossList)] <- paste0(i, '_', stepText)
    if(isTRUE(verbose)){
        print(paste0(date(), ': Iteration: ', i, ' ',
                stepText, ' loss: ', currLoss))
    }
    return(lossList)
}

lossED <- function(fds, lambda){
    K <- K(fds)
    N <- N(fds)
    rho <- matrix(rho(fds), ncol=ncol(K), nrow = nrow(K))
    D <- D(fds)

    y <- predictY(fds, noiseAlpha=currentNoiseAlpha(fds))

    return(fullNLL(y, rho, K, N, D, lambda))
}

# lossED <- function(fds){
#   K <- as.matrix(K(fds))
#   N <- as.matrix(N(fds))
#   mu <- predictMu(fds)
#   rho <- matrix(rho(fds), ncol=ncol(fds), nrow = nrow(fds))
#
#   r  <- (1-rho)/rho
#   eps <- 0.5
#   alpha  <- lgamma(mu*r)
#   alphaK <- lgamma(mu*r + K + eps)
#   beta   <- lgamma((mu-1)*(-r))
#   betaNK <- lgamma((mu-1)*(-r) + (N - K + eps))
#
#   #mean negative log likelihood with pseudocounts
#   mean(alpha + beta - alphaK - betaNK - lgamma(N+1+2*eps) + lgamma(K+1+eps) + lgamma(N-K+1+eps) + lgamma(r + N + 2*eps) - lgamma(r))
# }
