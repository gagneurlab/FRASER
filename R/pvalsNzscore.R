

# Function to calculate the z-scores
calculateZscore <- function(fds, type=currentType(fds)){
    currentType(fds) <- type
    
    counts(fds, type=type, side="other", HDF5=FALSE)      <- as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(counts(fds, type=type, side="ofInterest"))

    mu <- as.matrix(predictedMeans(fds))
    psi <- t(plogis(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE)))

    log2fc <- log2(psi) - log2(mu)

    # z = ( x - mean ) / sd
    zscores <- (log2fc - rowMeans(log2fc)) / rowSds(log2fc)

    zScores(fds) <- zscores

    return(fds)
}

# Function to calculate the p-values (both beta binomial and binomial)
calculatePvalues <- function(fds, type=currentType(fds),  BPPARAM=parallel(fds)){
    currentType(fds) <- type
    index <- getSiteIndex(fds, type=type)

    mu <- as.matrix(predictedMeans(fds))
    rho <- rho(fds)
    k <- as.matrix(K(fds))
    n <- as.matrix(N(fds))

    # beta binomial p-values
    pval_list <- bplapply(seq_len(nrow(mu)), singlePvalueBetaBinomial,
            k=k, n=n, mu=mu, rho=rho, BPPARAM=BPPARAM)
    pval <- do.call(rbind, pval_list)
    dval <- dbetabinom(k, n, mu, rho)
    pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
    fwer_pval <- bplapply(seq_len(ncol(pvals)), adjust_FWER_PValues, pvals=pvals, index, BPPARAM=BPPARAM)
    fwer_pvals <- do.call(cbind, fwer_pval)
    pVals(fds, dist="BetaBinom") <- fwer_pvals

    # binomial p-values
    pval_list <- bplapply(seq_len(nrow(mu)), singlePvalueBinomial, k=k, n=n, mu=mu, BPPARAM=BPPARAM)
    pval <- do.call(rbind, pval_list)
    dval <- dbinom(k, n , mu)
    pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
    fwer_pval <- bplapply(seq_len(ncol(pvals)), adjust_FWER_PValues, pvals=pvals, index, BPPARAM=BPPARAM)
    fwer_pvals <- do.call(cbind, fwer_pval)
    pVals(fds, dist="Binom") <- fwer_pvals

    return(fds)
}

adjust_FWER_PValues <- function(i, pvals=pvals, index=index){
    dt <- data.table(p=pvals[,i], idx=index)
    dt2 <- dt[,.(pa=min(p.adjust(p, method="holm"))),by=idx]
    setkey(dt2, "idx")[J(index)][,pa]
}

singlePvalueBetaBinomial <- function(idx, k, n, mu, rho){

    ki <- k[idx,]
    ni <- n[idx,]
    mui <- mu[idx,]
    rhoi <- rho[idx]

    pvals <- pmin(1, pbetabinom(ki, ni, mui, rhoi))
    return (pvals)
}

singlePvalueBinomial <- function(idx, k, n, mu){

  ki <- k[idx,]
  ni <- n[idx,]
  mui <- mu[idx,]

  pvals <- pmin(1, pbinom(ki, ni, mui))
  return (pvals)
}

calculatePadjValues <- function(fds, type=currentType(fds), method="BY"){
    currentType(fds) <- type
    index <- getSiteIndex(fds, type=type)
    idx   <- !duplicated(index)

    pvals <- pVals(fds, dist="BetaBinom")
    padj <- apply(pvals[idx,], 2, p.adjust, method=method)
    padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
    padjDT[,i:=NULL]
    padjVals(fds, dist="BetaBinom") <- as.matrix(padjDT)

    pvals <- pVals(fds, dist="Binom")
    padj <- apply(pvals[idx,], 2, p.adjust, method=method)
    padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
    padjDT[,i:=NULL]
    padjVals(fds, dist="Binom") <- as.matrix(padjDT)

    return(fds)
}

getSiteIndex <- function(fds, type){
    ans <- switch(type,
        psi5 = mcols(fds, type=type)[['startID']],
        psi3 = mcols(fds, type=type)[['endID']],
        psiSite = mcols(fds, type=type)[['spliceSiteID']]
    )
    return(ans)
}
