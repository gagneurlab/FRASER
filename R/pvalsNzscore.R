

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

# Function to calculate the p-values
calculatePvalues <- function(fds, type=currentType(fds), BPPARAM=bpparam()){
    currentType(fds) <- type

    mu <- as.matrix(predictedMeans(fds))
    rho <- matrix(rho(fds), nrow=nrow(mu) , ncol=ncol(mu))
    k <- as.matrix(K(fds))
    n <- as.matrix(N(fds))

    pval_list <- bplapply(seq_len(nrow(mu)), singlePvalue, k=k, n=n, mu=mu, rho=rho, BPPARAM=BPPARAM)
    pval <- matrix(unlist(pval_list), nrow=length(pval_list), byrow=TRUE)
    dval <- dbetabinom(k, n, mu, rho)
    pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)

    pVals(fds) <- pvals

    return(fds)
}

singlePvalue <- function(idx, k, n, mu, rho){

    ki <- k[idx,]
    ni <- n[idx,]
    mui <- mu[idx,]
    rhoi <- rho[idx]

    pvals <- pbetabinom(ki, ni, mui, rhoi)
    return (pvals)
}
