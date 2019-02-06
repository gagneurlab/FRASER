

# Function to calculate the z-scores
calculateZscore <- function(fds){

  mu <- predictedMeans(fds)
  psi <- t(plogis(x(fds)))
  log2fc <- log2(psi) - log2(mu)

  # z = ( x - mean ) / sd
  zscores <- (log2fc - rowMeans(log2fc)) / rowSds(log2fc)

  zScores(fds) <- zscores

  return(fds)

}

# Function to calculate the p-values
calculatePvalues <- function(fds){

  mu <- as.matrix(predictedMeans(fds))
  rho <- matrix(rho(fds), nrow=nrow(mu) , ncol=ncol(mu))
  k <- as.matrix(K(fds))
  n <- as.matrix(N(fds))

  pval <- pbetabinom(k, n, mu, rho)
  dval <- dbetabinom(k, n, mu, rho)
  pvals <- 2 * pmin(0.5, pval, 1 - pval + dval)

  pvalMat <- matrix(pvals, nrow=dim(mu)[1], ncol=dim(mu)[2])

  pVals(fds) <- pvalMat

  return(fds)
}
