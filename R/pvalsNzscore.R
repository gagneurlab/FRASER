

# Function to calculate the z-scores
calculateZscore <- function(fds, type=currentType(fds)){
    currentType(fds) <- type

    mu <- predictedMeans(fds)
    psi <- t(plogis(x(fds, all=TRUE)))

    log2fc <- log2(psi) - log2(mu)

    # z = ( x - mean ) / sd
    zscores <- (log2fc - rowMeans(log2fc)) / rowSds(log2fc)

    zScores(fds) <- zscores

    return(fds)
}

# Function to calculate the p-values (both beta binomial and binomial)
calculatePvalues <- function(fds, type=currentType(fds),  BPPARAM=bpparam()){
    currentType(fds) <- type

    mu <- as.matrix(predictedMeans(fds))
    rho <- rho(fds)
    k <- as.matrix(K(fds) +   pseudocount())
    n <- as.matrix(N(fds) + 2*pseudocount())

    # beta binomial p-values
    pval_list <- bplapply(seq_len(nrow(mu)), singlePvalueBetaBinomial, k=k, n=n, mu=mu, rho=rho, BPPARAM=BPPARAM)
    pval <- matrix(unlist(pval_list), nrow=length(pval_list), byrow=TRUE)
    dval <- dbetabinom(k, n, mu, rho)
    pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
    pVals(fds, dist="BetaBinom") <- pvals

    # binomial p-values
    pval_list <- bplapply(seq_len(nrow(mu)), singlePvalueBinomial, k=k, n=n, mu=mu, BPPARAM=BPPARAM)
    pval <- matrix(unlist(pval_list), nrow=length(pval_list), byrow=TRUE)
    dval <- dbinom(k, n , mu)
    pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
    pVals(fds, dist="Binom") <- pvals

    return(fds)
}

singlePvalueBetaBinomial <- function(idx, k, n, mu, rho){

    ki <- k[idx,]
    ni <- n[idx,]
    mui <- mu[idx,]
    rhoi <- rho[idx]

    pvals <- pbetabinom(ki, ni, mui, rhoi)
    return (pvals)
}

singlePvalueBinomial <- function(idx, k, n, mu){

  ki <- k[idx,]
  ni <- n[idx,]
  mui <- mu[idx,]

  pvals <- pmin(1, pbinom(ki, ni, mui))
  return (pvals)
}

calculatePadjValues <- function(fds, type=currentType(fds), method="BH"){
    currentType(fds) <- type

    pvals <- pVals(fds, dist="BetaBinom")
    padj <- apply(pvals, 2, p.adjust, method=method)
    padjVals(fds, dist="BetaBinom") <- padj

    pvals <- pVals(fds, dist="Binom")
    padj <- apply(pvals, 2, p.adjust, method=method)
    padjVals(fds, dist="Binom") <- padj

    return(fds)
}
