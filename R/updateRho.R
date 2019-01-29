#' 
#' Update theta step for autoencoder
#' 
#' @noRd
updateRho <- function(fds, rhoRange, BPPARAM, verbose){
  
  k <- as.matrix(counts(fds, type="psi3", side="ofInterest"))
  o <- as.matrix(counts(fds, type="psi3", side="other"))
  n <- k + o
  y <- predictY(fds)
  
  #fitparameters <- bplapply(seq_along(fds), estRho, 
  #                          k=k, n=n, mu=mu, rhoRange=rhoRange,
  #                          BPPARAM=BPPARAM, nll=negLogLikelihoodRho)
  
  fitparameters <- bplapply(seq_along(fds), estRho, 
                            k=k, n=n, y=y, rhoRange=rhoRange,
                            BPPARAM=BPPARAM, nll=truncNLL_rho)
  
  rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
  
  if(isTRUE(verbose)){
    print(summary(rho(fds)))
  }
  
  validObject(fds)
  return(fds)
}

estRho <- function(idx, k, n, y, rhoRange, nll, control=list()){
  ki <- k[idx,]
  ni <- n[idx,]
  yi <- y[idx,]
  
  est <- optimize(f=nll, interval=rhoRange, yi=yi, ki=ki, ni=ni, maximum = FALSE, tol=0.0000001)
}

negLogLikelihoodRho <- function(rho, ki, ni, mui){
  #-mean(dbetabinom(ki + 0.5, ni + 1, mu, rho, log=TRUE))
  
  r  <- (1-rho)/rho
  eps <- 0.5
  alpha  <- lgamma(mui*r) 
  alphaK <- lgamma(mui*r + ki + eps) 
  beta   <- lgamma((mui-1)*(-r)) 
  betaNK <- lgamma((mui-1)*(-r) + (ni - ki + eps)) 
  
  #mean negative log likelihood with pseudocounts
  mean(alpha + beta - alphaK - betaNK - lgamma(ni+1+2*eps) + lgamma(ki+1+eps) + lgamma(ni-ki+1+eps) + lgamma(r + ni + 2*eps) - lgamma(r))
}

trunc_negLogLikelihoodRho <- function(rho, ki, ni, mui){
  #-mean(dbetabinom(ki + 0.5, ni + 1, mu, rho, log=TRUE))
  
  r  <- (1-rho)/rho
  eps <- 0.5
  alpha  <- lgamma(mui*r) 
  alphaK <- lgamma(mui*r + ki + eps) 
  beta   <- lgamma((mui-1)*(-r)) 
  betaNK <- lgamma((mui-1)*(-r) + (ni - ki + eps)) 
  
  #mean negative log likelihood with pseudocounts
  mean(alpha + beta - alphaK - betaNK )
}
