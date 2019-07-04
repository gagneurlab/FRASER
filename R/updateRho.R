#'
#' Update theta step for autoencoder
#'
#' @noRd
updateRho <- function(fds, type, rhoRange, BPPARAM, verbose){

  k <- K(fds)
  n <- N(fds)
  y <- predictY(fds, noiseAlpha=currentNoiseAlpha(fds))

  fitparameters <- bplapply(seq_len(nrow(k)), estRho,
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


methodOfMomentsRho <- function(k, n, rhoRange=c(1e-5, 1 - 1e-5)){
  # taken from wiki: https://en.wikipedia.org/wiki/Beta-binomial_distribution#Method_of_moments
  k <- as.matrix(k) + pseudocount()
  mode(k) <- "double"
  n <- as.matrix(n) + 2*pseudocount()
  mode(n) <- "double"
  
  # normalizing counts with respect to mean(n) because for each junction an equal n for all samples is needed
  nM <- round(rowMeans(n))
  kNorm <- round(nM/n * k)
  
  N  <- ncol(kNorm)
  m1 <- rowSums(kNorm)/N
  m2 <- rowSums(kNorm*kNorm)/N
  
  denom <- (nM*(m2/m1 - m1 - 1) + m1)
  a <- (nM*m1 - m2) / denom
  b <- ((nM - m1)*(nM - m2/m1))/ denom
  rho <- 1 / (1 + a + b)
  
  # set correct boundaries
  rho[a < 0 | b < 0] <- rhoRange[1]
  rho[is.na(rho)] <- rhoRange[1]
  rho <- pmax(rho, rhoRange[1])
  rho <- pmin(rho, rhoRange[2])
  
  return(rho)
}
