
#' @describeIn FRASER This function calculates z-scores based on the 
#' observed and expected logit 
#' psi.
#' 
#' @param logit Indicates if z scores are computed on the logit scale (default) 
#'      or in the natural (psi) scale.
#' @export
calculateZscore <- function(fds, type=currentType(fds), logit=TRUE){
    currentType(fds) <- type

    mu <- predictedMeans(fds)
    psi <- (K(fds) + pseudocount()) / (N(fds) + 2*pseudocount())
    
    if(isTRUE(logit)){
        mu <- qlogis(mu)
        psi <- qlogis(psi)
    }
    
    if(is.matrix(psi) && !is.matrix(mu)){
        mu <- as.matrix(mu)
    }

    residual <- psi - mu

    # z = ( x - mean ) / sd
    zscores <- (residual - rowMeans(residual)) / rowSds(residual)

    zScores(fds, withDimnames=FALSE) <- zscores
    
    return(fds)
}

#' @describeIn FRASER This function calculates two-sided p-values based on 
#' the beta-binomial distribution (or binomial or normal if desired). The 
#' returned p values are not yet adjusted with Holm's method per donor or 
#' acceptor site, respectively. 
#' 
#' @param distributions The distribution based on which the p-values are 
#' calculated. Possible are beta-binomial, binomial and normal.
#' @param capN Counts are capped at this value to speed up the p-value 
#' calculation
#' 
#' @export
calculatePvalues <- function(fds, type=currentType(fds),
                    implementation="PCA", BPPARAM=bpparam(),
                    distributions=c("betabinomial"), capN=5*1e5){
    distributions <- match.arg(distributions, several.ok=TRUE,
            choices=c("betabinomial", "binomial", "normal"))
    
    # make sure its only in-memory data for k and n
    currentType(fds) <- type
    fds <- putCounts2Memory(fds, type)
    
    # if method BB is used take the old FRASER code
    if(implementation %in% c("BB")){
        index <- getSiteIndex(fds, type)
        pvals <- getAssayMatrix(fds, "pvalues_BB", type)
        fwer_pval  <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                pvals=pvals, index, BPPARAM=BPPARAM)
        fwer_pvals <- do.call(cbind, fwer_pval)
        pVals(fds, type=type, dist="BetaBinomial", 
                withDimnames=FALSE) <- fwer_pvals
        return(fds)
    }
    
    index <- getSiteIndex(fds, type=type)
    mu <- as.matrix(predictedMeans(fds))
    rho <- rho(fds)
    alpha <- mu * (1 - rho)/rho
    beta <- (1 - mu) * (1 - rho)/rho
    k <- K(fds)
    n <- N(fds)
    
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
        # # above code would be nicer but fails for assignment to 
        # # delayedMatrix, for which the following is needed:
        # bigN <- which(n > capN, arr.ind=TRUE)
        # if(nrow(bigN) >= 1){
        #     for(ind in seq_len(nrow(bigN))){
        #         i <- bigN[ind, 1]
        #         j <- bigN[ind, 2]
        #         facN <- capN/n[i,j]
        #         k[i,j] <- pmin(round(k[i,j] * facN), capN)
        #         n[i,j] <- capN
        #     }
        # }
    }
    
    if("betabinomial" %in% distributions){
        # beta binomial p-values
        pval_list <- bplapply(seq_row(mu), singlePvalueBetaBinomial,
                k=k, n=n, mu=mu, rho=rho, BPPARAM=BPPARAM)
        pval <- do.call(rbind, pval_list)
        dval <- matrix(nrow=nrow(k), ncol=ncol(k), 
                        dbbinom(as.matrix(k), as.matrix(n), 
                                as.matrix(alpha), as.matrix(beta)))
        pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
        pVals(fds, dist="BetaBinomial", level="junction", 
                withDimnames=FALSE) <- pvals
    }
    
    if("binomial" %in% distributions){
        # binomial p-values
        pval_list <- bplapply(seq_row(mu), singlePvalueBinomial, k=k, n=n,
                mu=mu, BPPARAM=BPPARAM)
        pval <- do.call(rbind, pval_list)
        dval <- dbinom(as.matrix(k), as.matrix(n), as.matrix(mu))
        pvals <- 2 * pmin(pval, 1 - pval + dval, 0.5)
        pVals(fds, dist="Binomial", level="junction", 
                withDimnames=FALSE) <- pvals
    }
    
    if("normal" %in% distributions){
        fds <- putCounts2Memory(fds, type)
        yin <- t(x(fds, all=TRUE, noiseAlpha=NULL, center=FALSE))
        yout <- predictY(fds, type, noiseAlpha=NULL)
        epsilon <- yin - yout
        rsd <- rowSds(epsilon)
        pval <- pnorm(epsilon, sd=rsd)
        pvals <- 2 * pmin(pval, 1 - pval, 0.5)
        pVals(fds, dist="Normal", level="junction", 
                withDimnames=FALSE) <- pvals
    }
    
    fds
}

adjust_FWER_PValues <- function(i, pvals, index, rho, rhoCutoff, 
                                method="holm"){
    dt <- data.table(p=pvals[,i], idx=index, rho=rho)
    dt[rho > rhoCutoff, p:=NA]
    suppressWarnings(dt2 <- dt[,.(pa=min(p.adjust(p, method=method), 
                                        na.rm=TRUE)),by=idx])
    dt2[is.infinite(pa), pa:=NA]
    setkey(dt2, "idx")[J(index)][,pa]
}

adjust_FWER_PValues_per_idx <- function(i, pvals, index, rho, rhoCutoff, 
                                method="holm"){
    pvals[rho > rhoCutoff,] <- NA
    dttmp <- data.table(idx=index, rho=rho, 
                        apply(pvals, 2, as.numeric))[idx == i,]
    suppressWarnings(
        pa <- apply(as.matrix(dttmp[,-c("idx", "rho")]), 2,
                        function(x) min(p.adjust(x, method=method), 
                                        na.rm = TRUE) )
    )
    pa[is.infinite(pa)] <- NA
    return(pa)
}

getFWERpvals_bySample <- function(pvals, index, rho, method="holm", 
                                rhoCutoff=0.1, BPPARAM=bpparam()){
    fwer_pval <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                    pvals=pvals, index, BPPARAM=BPPARAM,
                    method=method, rho=rho, rhoCutoff=rhoCutoff)
    fwer_pvals <- do.call(cbind, fwer_pval)
    return(fwer_pvals)
}

getFWERpvals_byIdx <- function(pvals, index, rho, method="holm", 
                                rhoCutoff=0.1, BPPARAM=bpparam()){
    unique_idx <- unique(index)
    fwer_pval <- bplapply(unique_idx, adjust_FWER_PValues_per_idx,
                    pvals=pvals, index, BPPARAM=BPPARAM,
                    method=method, rho=rho, rhoCutoff=rhoCutoff)
    fwer_pvals <- do.call(rbind, fwer_pval)
    fwer_pvals <- as.matrix(
        setkey(data.table(idx=unique_idx, fwer_pvals), 
                "idx")[J(index)][,-c("idx")])
    return(fwer_pvals)
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
        message(date(), ": obtained NA pvalues for junction ", idx)
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

#' @describeIn FRASER This function adjusts the previously calculated 
#' p-values per sample for multiple testing. First, the previoulsy calculated 
#' junction-level p values are adjusted with Holm's method per donor or 
#' acceptor site, respectively. Then, if gene symbols have been annotated to 
#' junctions (and not otherwise requested), gene-level p values are computed.
#' 
#' @param method The p.adjust method that should be used for genome-wide 
#' multiple testing correction.
#' @param rhoCutoff The cutoff value on the fitted rho value 
#'     (overdispersion parameter of the betabinomial) above which junctions are 
#'     masked with NA during p value adjustment. 
#' @param geneLevel Logical value indiciating whether gene-level p values 
#'     should be calculated. Defaults to TRUE.
#' @param geneColumn The column name of the column that has the gene annotation 
#'              that will be used for gene-level pvalue computation.
#' 
#' @export
calculatePadjValues <- function(fds, type=currentType(fds), method="BY",
                                rhoCutoff=0.1, geneLevel=TRUE, 
                                geneColumn="hgnc_symbol", BPPARAM=bpparam()){
    currentType(fds) <- type
    index <- getSiteIndex(fds, type=type)
    idx   <- !duplicated(index)
    
    for(i in c("BetaBinomial", "Binomial", "Normal")){
        # only do it if it exists
        if(!paste0("pvalues", i, "_junction_", type) %in% assayNames(fds)){
            next
        }
        
        pvals <- pVals(fds, dist=i, level="junction")
        rho <- rho(fds, type=type)
        
        # splice site-level pval correction
        message(date(), ":  adjusting junction-level pvalues ...")
        # fwer_pvals <- getFWERpvals_byIdx(pvals, index, rho, method="holm", 
        #                                 rhoCutoff=rhoCutoff, BPPARAM=BPPARAM)
        fwer_pvals <- getFWERpvals_bySample(pvals, index, rho, method="holm", 
                                        rhoCutoff=rhoCutoff, BPPARAM=BPPARAM)
        pVals(fds, dist=i, level="site", filters=list(rho=rhoCutoff),
                withDimnames=FALSE) <- fwer_pvals
        
        # junction-level FDR correction
        message(date(), ":  genome-wide FDR for junction-level pvalues ...")
        padj <- apply(fwer_pvals[idx,], 2, p.adjust, method=method)
        padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
        padjDT[,i:=NULL]
        padjVals(fds, dist=i, level="site", filters=list(rho=rhoCutoff),
                withDimnames=FALSE) <- as.matrix(padjDT)
        
        # gene-level pval correction and FDR
        if(geneColumn %in% colnames(mcols(fds, type=type)) && 
                isTRUE(geneLevel)){
            message(date(), ":  calculating gene-level pvalues ...")
            gene_pvals <- getPvalsPerGene(fds=fds, type=type, pvals=fwer_pvals,
                                            method="holm", FDRmethod=method, 
                                            geneColumn=geneColumn,
                                            BPPARAM=BPPARAM)
            pVals(fds, dist=i, level="gene", filters=list(rho=rhoCutoff),
                    withDimnames=FALSE) <- gene_pvals[["pvals"]]
            padjVals(fds, dist=i, level="gene", filters=list(rho=rhoCutoff),
                    withDimnames=FALSE) <- gene_pvals[["padj"]]
        } else{
            warning("Gene-level pvalues could not be calculated as column ",
                    geneColumn, "\nwas not found for the given fds object. ", 
                    "Please annotate gene symbols \nfirst using the ", 
                    "annotateRanges function.")
        }
    }
    
    return(fds)
}

getPvalsPerGene <- function(fds, type, 
                    pvals=pVals(fds, type=type, level="site"),
                    sampleID=NULL, method="holm", FDRmethod="BY", 
                    geneColumn="hgnc_symbol", BPPARAM=bpparam()){
    # extract data and take only the first index of per site
    message(date(), ":   starting gene-level pval computation for type ", type)
    samples <- samples(fds)
    if(is.null(colnames(pvals))){
        colnames(pvals) <- samples
    }
    dt <- data.table(
            idx=getSiteIndex(fds, type=type),
            geneID=getGeneIDs(fds, type=type, unique=FALSE, 
                                geneColumn=geneColumn),
            as.data.table(pvals))
    dt <- dt[!is.na(geneID)]
    geneIDs <- getGeneIDs(fds, type=type, unique=TRUE, 
                            geneColumn=geneColumn)
    
    # separate geneIDs into rows
    dt[, dt_idx:=seq_len(.N)]
    dt_tmp <- dt[, splitGenes(geneID), by="dt_idx"]
    dt <- dt[dt_tmp$dt_idx,]
    dt[,`:=`(geneID=dt_tmp$V1, dt_idx=NULL)]
    setkey(dt, geneID)
    
    # extract samples
    if(!is.null(sampleID)){
        samples <- sampleID
    }
    
    # aggregate pvalues to gene level per sample
    message(date(), ":   gene-level pval computation per gene (n=", 
            length(geneIDs), ")")
    pvalsPerGene <- genePvalsByGeneID(dt, samples=samples, geneIDs=geneIDs, 
                                        method=method, BPPARAM=BPPARAM)
    
    # compute FDR
    message(date(), ":   genome-wide FDR for gene-level pvals for type ", type)
    padjPerGene <- apply(pvalsPerGene, 2, p.adjust, method=FDRmethod)
    
    message(date(), ":   finished gene-level pval computation for type ", type)
    return(list(pvals=pvalsPerGene, padj=padjPerGene))

}

getSiteIndex <- function(fds, type){
    if(type == "theta"){
        return(mcols(fds, type=type)[['spliceSiteID']])
    }
    
    if(type == "jaccard"){
        return(seq_len(nrow(fds)))
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

getGeneIDs <- function(fds, type, unique=TRUE, geneColumn="hgnc_symbol"){
    geneIDs <- mcols(fds, type=type)[[geneColumn]]
    if(isTRUE(unique)){
        geneIDs <- unique(unlist(lapply(geneIDs, FUN=function(g){
            unlist(strsplit(g, ";"))}) ))
        geneIDs <- geneIDs[!is.na(geneIDs)]
    }
    geneIDs
}

genePvalsByGeneID <- function(dt, samples, geneIDs, method, BPPARAM){
    pvalsPerGene <- bplapply(geneIDs, BPPARAM=BPPARAM,
        FUN=function(g) {
            dttmp <- dt[geneID == g][!duplicated(idx)]
            suppressWarnings(
                pval_g <- apply(as.matrix(dttmp[,-c("idx", "geneID")]), 2,
                    function(x) min(p.adjust(x, method=method), na.rm = TRUE) )
            )
            pval_g[is.infinite(pval_g)] <- NA
            pval_g
        })
    pvalsPerGene <- do.call(rbind, pvalsPerGene)
    rownames(pvalsPerGene) <- geneIDs
    return(pvalsPerGene)
}

