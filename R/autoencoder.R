#' 
#' Main autoencoder fit function
#' 
#' @noRd
fitAutoencoder <- function(fds, q, rhoRange=c(1e-6, 1-1e-6), lambda=10,
                           convergence=1e-5, iterations=15, initialize=TRUE,
                           control=list(), BPPARAM=bpparam(), verbose=FALSE){
  
  # Check input

  
  if(!bpisup(BPPARAM)){
    bpstart(BPPARAM)
  }
  
  # initialize E and D using PCA and bias as zeros.
  if(isTRUE(initialize) | is.null(E(fds)) | is.null(D(fds))){
    fds <- initAutoencoder(fds, q, rhoRange)
  }
  
  # initial loss
  lossList <- c(init_pca=lossED(fds, lambda))
  print(paste0('Initial PCA loss: ', lossList[1]))
  
  #inital rho values
  if(isTRUE(verbose)){
    print(summary(rho(fds)))
  }
  
  # initialize D 
  fds <- updateD(fds, lambda=lambda, control=control, BPPARAM=BPPARAM, verbose=verbose)
  lossList <- updateLossList(fds, lossList, 'init', 'D', lambda, verbose=verbose)
  
  # initialize rho step
  fds <- updateRho(fds, rhoRange, BPPARAM=BPPARAM, verbose=verbose)
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
    fds <- updateD(fds, lambda=lambda, control=control, BPPARAM=BPPARAM, verbose=verbose)
    lossList <- updateLossList(fds, lossList, i, 'D', lambda, verbose=verbose)
    
    # update rho step
    fds <- updateRho(fds, rhoRange, BPPARAM=BPPARAM, verbose=verbose)
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
  
  bpstop(BPPARAM)
  print(Sys.time() - t1)
  
  print(paste0(i, ' Final betabin-AE loss: ', lossList[length(lossList)]))
  
  
  # add additional values for the user to the object
  metadata(fds)[['dim']] <- dim(fds)
  metadata(fds)[['loss']] <- lossList
  metadata(fds)[['convList']] <- lossList
  
  validObject(fds)
  return(fds)
}

initAutoencoder <- function(fds, q, rhoRange){
  
  pca <- pca(x(fds), nPcs=q)
  pc  <- loadings(pca)
  
  # Set initial values from PCA
  D(fds) <- pc
  E(fds) <- pc
  b(fds) <- double(nrow(fds))
  
  # initialize rho
  rho(fds) <- rowSds(plogis(t(x(fds))))
  rho(fds)[rho(fds) < rhoRange[1]] <- rhoRange[1]
  
  # reset counters 
  mcols(fds)['NumConvergedD'] <- 0
  
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
  K <- as.matrix(K(fds))
  N <- as.matrix(N(fds))
  mu <- predictMu(fds)
  rho <- matrix(rho(fds), ncol=ncol(fds), nrow = nrow(fds))
  #rho <- matrix(rho(fds), nrow=ncol(fds), ncol = nrow(fds), byrow = TRUE)
  D <- D(fds)

  return(fullNLL(mu, rho, K, N, D, lambda))
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
