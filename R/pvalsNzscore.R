

# Function to calculate the z-scores
calculateZscore <- function(fds){
  
  mu <- predictMu(fds, all=TRUE)
  psi <- t(plogis(x(fds, all=TRUE)))
  log2fc <- log2(psi) - log2(mu)
  
  # z = ( x - mean ) / sd
  zscores <- (log2fc - rowMeans(log2fc)) / rowSds(log2fc)
  
  zScores(fds) <- zscores
  
  return(fds)
  
}

# Function to calculate the p-values
calculatePvalues <- function(fds){
  
  mu <- predictMu(fds, all=TRUE)
  rho <- matrix(rho(fds), nrow = nrow(mu) , ncol = ncol(mu))
  k <- as.matrix(K(fds))
  n <- as.matrix(N(fds))
  
  pval <- VGAM::pbetabinom(k, n, mu, rho)
  dval <- VGAM::dbetabinom(k, n, mu, rho)
  pvals <- 2 * pmin(0.5, pval + dval, 1 - pval + dval)
  
  pvalMat <- matrix(pvals, nrow=dim(fds)[1], ncol=dim(fds)[2])
  
  pVals(fds) <- pvalMat
  
  return(fds)
}
