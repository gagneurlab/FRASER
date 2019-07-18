

# Function to calculate the z-scores
calculateZscore <- function(fds, type=currentType(fds), correction="FraseR"){
    currentType(fds) <- type

    counts(fds, type=type, side="other", HDF5=FALSE)      <-
            as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <-
            as.matrix(counts(fds, type=type, side="ofInterest"))

    logit_mu <- qlogis(as.matrix(predictedMeans(fds)))
    logit_psi <- t(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE))

    logitfc <- logit_psi - logit_mu

    # z = ( x - mean ) / sd
    zscores <- (logitfc - rowMeans(logitfc)) / rowSds(logitfc)

    zScores(fds) <- zscores

    return(fds)
}

# Function to calculate the p-values (both beta binomial and binomial)
calculatePvalues <- function(fds, type=currentType(fds),
                    correction="FraseR", BPPARAM=parallel(fds)){

    # make sure its only in-memory data for k and n
    currentType(fds) <- type
    counts(fds, type=type, side="other", HDF5=FALSE) <- as.matrix(N(fds) - K(fds))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(K(fds))

    # if method BB is used take the old FraseR code
    if(correction %in% c("BB")){
        index <- getSiteIndex(fds, type)
        pvals <- getAssayMatrix(fds, "pvalues_BB", type)
        fwer_pval  <- bplapply(seq_len(ncol(pvals)), adjust_FWER_PValues,
                pvals=pvals, index, BPPARAM=BPPARAM)
        fwer_pvals <- do.call(cbind, fwer_pval)
        pVals(fds, type=type, dist="BetaBinomial") <- fwer_pvals
        return(fds)
    }

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
    pVals(fds, dist="BetaBinomial") <- fwer_pvals

    # binomial p-values
    pval_list <- bplapply(seq_len(nrow(mu)), singlePvalueBinomial, k=k, n=n, mu=mu, BPPARAM=BPPARAM)
    pval <- do.call(rbind, pval_list)
    dval <- dbinom(k, n , mu)
    pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
    fwer_pval <- bplapply(seq_len(ncol(pvals)), adjust_FWER_PValues, pvals=pvals, index, BPPARAM=BPPARAM)
    fwer_pvals <- do.call(cbind, fwer_pval)
    pVals(fds, dist="Binomial") <- fwer_pvals

    return(fds)
}

adjust_FWER_PValues <- function(i, pvals=pvals, index=index){
    dt <- data.table(p=pvals[,i], idx=index)
    dt2 <- dt[,.(pa=min(p.adjust(p, method="holm"), na.rm=TRUE)),by=idx]
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

    pvals <- pVals(fds, dist="BetaBinomial")
    padj <- apply(pvals[idx,], 2, p.adjust, method=method)
    padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
    padjDT[,i:=NULL]
    padjVals(fds, dist="BetaBinomial") <- as.matrix(padjDT)

    # only do it if it exists
    if(paste0("pvaluesBinomial_", type) %in% assayNames(fds)){
        pvals <- pVals(fds, dist="Binomial")
        padj <- apply(pvals[idx,], 2, p.adjust, method=method)
        padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
        padjDT[,i:=NULL]
        padjVals(fds, dist="Binomial") <- as.matrix(padjDT)
    }

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

getGeneIDs <- function(fds, type){
    geneIDs <- unique(mcols(fds, type=type)$hgnc_symbol)
    geneIDs <- geneIDs[!is.na(geneIDs)]
    return(geneIDs)
}

getPvalsPerGene <- function(fds, type, pvals=pVals(fds, type=type), sampleID=NULL,
                            method="holm", BPPARAM=parallel(fds)){

  # extract data
  dt <- data.table(
    geneID=mcols(fds, type=type)$hgnc_symbol,
    as.data.table(pvals) )
  setkey(dt, geneID)

  samples <- samples(fds)
  if(!is.null(sampleID)){
    samples <- sampleID
  }

  pvalsPerGene <- matrix(unlist(bplapply(samples, BPPARAM=BPPARAM, function(i){
    dttmp <- dt[,min(p.adjust(get(i), method=method)),by=geneID]
    dttmp <- dttmp[!is.na(geneID),]
    setkey(dttmp, geneID)
    dttmp[J(getGeneIDs(fds, type)), V1]
  })), ncol=length(samples))

  return(pvalsPerGene)

}

getPadjPerGene <- function(pvals, method="BY"){

  padjPerGene <- apply(pvals, 2, p.adjust, method=method)

  return(padjPerGene)

}
