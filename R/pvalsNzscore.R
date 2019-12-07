
#' 
#' @title Calculate P-values and Z-scores
#'  
#' @description The FraseR package provides functions to calculate nominal and 
#' and adjusted p-values as well as z-scores after fitting the beta-binomial 
#' parameters. See Detail and Functions for more information.
#'  
#' @details See Functions.
#' 
#' @inheritParams fit 
#' 
#' @return FraseRDataSet
#' 
#' @name pvalZscore
#' @rdname pvalZscore
#' 
#' @examples
#'   # preprocessing
#'   fds <- createTestFraseRSettings()
#'   fds <- countRNAData(fds)
#'   fds <- calculatePSIValues(fds)
#'   fds <- fit(fds, q=5, correction="PCA", type="psi5")
#'   
#'   # nomial p values
#'   fds <- calculatePvalues(fds, type="psi5")
#'   head(pVals(fds, type="psi5"))
#'   
#'   # donor site adjusted p values
#'   fds <- calculatePadjValues(fds, type="psi5", method="BY")
#'   head(padjVals(fds, type="psi5"))
#'   
#'   # z scores
#'   fds <- calculateZscore(fds, type="psi5")
#'   head(zScores(fds, type="psi5")) 
#'
NULL

#' @describeIn pvalZscore This function calculates z-scores based on the 
#' observed and expected logit 
#' psi.
#' 
#' @export
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

#' @describeIn pvalZscore This function calculates two-sided p-values based on 
#' the beta-binomial distribution (or binomial or normal if desired).
#' 
#' @param distributions The distribution based on which the p-values are 
#' calculated. Possible are beta-binomial, binomial and normal.
#' @param capN Counts are capped at this value to speed up the p-value 
#' calculation
#' 
#' @export
calculatePvalues <- function(fds, type=currentType(fds),
                    correction="FraseR", BPPARAM=bpparam(),
                    distributions=c("betabinomial"), capN=5*1e5){
    distributions <- match.arg(distributions, several.ok=TRUE,
            choices=c("betabinomial", "binomial", "normal"))
    
    # make sure its only in-memory data for k and n
    currentType(fds) <- type
    fds <- putCounts2Memory(fds, type)
    
    # if method BB is used take the old FraseR code
    if(correction %in% c("BB")){
        index <- getSiteIndex(fds, type)
        pvals <- getAssayMatrix(fds, "pvalues_BB", type)
        fwer_pval  <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                pvals=pvals, index, BPPARAM=BPPARAM)
        fwer_pvals <- do.call(cbind, fwer_pval)
        pVals(fds, type=type, dist="BetaBinomial") <- fwer_pvals
        return(fds)
    }
    
    index <- getSiteIndex(fds, type=type)
    mu <- as.matrix(predictedMeans(fds))
    rho <- rho(fds)
    alpha <- mu * (1 - rho)/rho
    beta <- (1 - mu) * (1 - rho)/rho
    k <- as.matrix(K(fds))
    n <- as.matrix(N(fds))
    
    # betaBinomial functions get slowed down drastically if
    # N is big (2 mio and bigger). Hence downsample if requested to 1mio max
    if(isTRUE(capN)){
        capN <- 1e6
    }
    if(isScalarNumeric(capN)){
        bigN <- which(n > capN)
        if(length(bigN) >= 1){
            facN <- capN/n[bigN]
            k[bigN] <- pmin(round(k[bigN] * facN), capN)
            n[bigN] <- capN
        }
    }
    
    if("betabinomial" %in% distributions){
        # beta binomial p-values
        pval_list <- bplapply(seq_row(mu), singlePvalueBetaBinomial,
                k=k, n=n, mu=mu, rho=rho, BPPARAM=BPPARAM)
        pval <- do.call(rbind, pval_list)
        dval <- matrix(dbbinom(k, n, alpha, beta), nrow=nrow(k), ncol=ncol(k))
        pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
        fwer_pval <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                pvals=pvals, index, BPPARAM=BPPARAM)
        fwer_pvals <- do.call(cbind, fwer_pval)
        pVals(fds, dist="BetaBinomial") <- fwer_pvals
    }
    
    if("binomial" %in% distributions){
        # binomial p-values
        pval_list <- bplapply(seq_row(mu), singlePvalueBinomial, k=k, n=n,
                mu=mu, BPPARAM=BPPARAM)
        pval <- do.call(rbind, pval_list)
        dval <- dbinom(k, n , mu)
        pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
        fwer_pval <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                pvals=pvals, index, BPPARAM=BPPARAM)
        fwer_pvals <- do.call(cbind, fwer_pval)
        pVals(fds, dist="Binomial") <- fwer_pvals
    }
    
    if("normal" %in% distributions){
        yin <- t(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE))
        yout <- predictY(fds, type, noiseAlpha=NULL)
        epsilon <- yin - yout
        rsd <- rowSds(epsilon)
        pval <- pnorm(epsilon, sd=rsd)
        pvals <- 2 * pmin(pval, 1 - pval, 0.5)
        fwer_pval <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                pvals=pvals, index, BPPARAM=BPPARAM)
        fwer_pvals <- do.call(cbind, fwer_pval)
        pVals(fds, dist="Normal") <- fwer_pvals
    }
    
    fds
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
    alphai <- mui * (1 - rhoi)/rhoi
    betai <- (1 - mui) * (1 - rhoi)/rhoi

    # try catch block to overcome long running times
    pvals <- pmin(1, pbbinom(ki, ni, alphai, betai))

    if(any(is.na(pvals))){
        message(date(), " : ", idx)
    }

    return (pvals)
}

singlePvalueBinomial <- function(idx, k, n, mu){

    ki <- k[idx,]
    ni <- n[idx,]
    mui <- mu[idx,]

    pvals <- pmin(1, pbinom(ki, ni, mui))
    return (pvals)
}

#' @describeIn pvalZscore This function adjusts the previously calculated 
#' p-values per donor/acceptor site for multiple testing.
#' 
#' @param method The p.adjust method that should be used. 
#' 
#' @export
calculatePadjValues <- function(fds, type=currentType(fds), method="BY"){
    currentType(fds) <- type
    index <- getSiteIndex(fds, type=type)
    idx   <- !duplicated(index)
    
    for(i in c("BetaBinomial", "Binomial", "Normal")){
        # only do it if it exists
        if(!paste0("pvalues", i, "_", type) %in% assayNames(fds)){
            next
        }
        
        pvals <- pVals(fds, dist=i)
        padj <- apply(pvals[idx,], 2, p.adjust, method=method)
        padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
        padjDT[,i:=NULL]
        padjVals(fds, dist=i) <- as.matrix(padjDT)
    }
    
    return(fds)
}

getSiteIndex <- function(fds, type){
    if(type == "psiSite"){
        return(mcols(fds, type=type)[['spliceSiteID']])
    }
    
    startId <- mcols(fds, type=type)[,"startID"]
    endId   <- mcols(fds, type=type)[,"endID"]
    strand  <- strand(rowRanges(fds, type=type))
    strand[strand == "*"] <- "+"
    
    selectionMat <- as.matrix(data.frame(row=seq_along(startId), 
            col=1 + as.vector(
                    type == "psi5" & strand == "-" | 
                    type == "psi3" & strand == "+")))
    
    ans <- as.matrix(cbind(startId, endId))
    ans[selectionMat]
}

getGeneIDs <- function(fds, type, unique=TRUE){
    geneIDs <- mcols(fds, type=type)$hgnc_symbol
    if(isTRUE(unique)){
        geneIDs <- unique(geneIDs)
        geneIDs <- geneIDs[!is.na(geneIDs)]
    }
    geneIDs
}

getPvalsPerGene <- function(fds, type, pvals=pVals(fds, type=type),
                    sampleID=NULL, method="holm", BPPARAM=bpparam()){
    
    # extract data and take only the first index of per site
    dt <- data.table(
            index=getSiteIndex(fds, type=type),
            geneID=getGeneIDs(fds, type=type, unique=FALSE),
            as.data.table(pvals))[!duplicated(index)]
    dt <- dt[!is.na(geneID)]
    setkey(dt, geneID)
    
    samples <- samples(fds)
    if(!is.null(sampleID)){
        samples <- sampleID
    }
    
    pvalsPerGene <- matrix(unlist(bplapply(samples, BPPARAM=BPPARAM,
        function(i){
                dttmp <- dt[,min(p.adjust(get(i), method=method)),by=geneID]
                setkey(dttmp, geneID)
                dttmp[J(getGeneIDs(fds, type=type)), V1]
        })), ncol=length(samples))
    
    colnames(pvalsPerGene) <- samples
    rownames(pvalsPerGene) <- getGeneIDs(fds, type=type)
    
    return(pvalsPerGene)

}

getPadjPerGene <- function(pvals, method="BY"){

    padjPerGene <- apply(pvals, 2, p.adjust, method=method)

    return(padjPerGene)

}
