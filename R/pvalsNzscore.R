
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
                                rhoCutoff, BPPARAM=bpparam()){
    fwer_pval <- bplapply(seq_col(pvals), adjust_FWER_PValues,
                    pvals=pvals, index, BPPARAM=BPPARAM,
                    method=method, rho=rho, rhoCutoff=rhoCutoff)
    fwer_pvals <- do.call(cbind, fwer_pval)
    return(fwer_pvals)
}

getFWERpvals_byIdx <- function(pvals, index, rho, method="holm", 
                                rhoCutoff, BPPARAM=bpparam()){
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
#'     masked with NA during p value adjustment (default: NA, no masking). 
#' @param geneLevel Logical value indiciating whether gene-level p values 
#'     should be calculated. Defaults to TRUE.
#' @param geneColumn The column name of the column that has the gene annotation 
#'              that will be used for gene-level pvalue computation.
#' @param subsets A named list of named lists specifying any number of gene 
#'             subsets (can differ per sample). For each subset, FDR correction
#'             will be limited to genes in the subset, and the FDR corrected 
#'             pvalues stored as an assay in the fds object in addition to the 
#'             transcriptome-wide FDR corrected pvalues. See the examples for 
#'             how to use this argument.
#' 
#' @export
calculatePadjValues <- function(fds, type=currentType(fds), method="BY",
                                rhoCutoff=NA, geneLevel=TRUE, 
                                geneColumn="hgnc_symbol", subsets=NULL, 
                                BPPARAM=bpparam()){
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
        fwer_pvals <- getFWERpvals_bySample(pvals, index, rho, method="holm", 
                            rhoCutoff=ifelse(is.na(rhoCutoff), 1, rhoCutoff), 
                            BPPARAM=BPPARAM)
        if(!is.na(rhoCutoff)){
            filters <- list(rho=rhoCutoff)
        } else{
            filters <- list()
        }
        pVals(fds, dist=i, level="site", filters=filters,
                withDimnames=FALSE) <- fwer_pvals
        
        # junction-level FDR correction
        message(date(), ":  genome-wide FDR for junction-level pvalues ...")
        padj <- apply(fwer_pvals[idx,], 2, p.adjust, method=method)
        padjDT <- data.table(cbind(i=unique(index), padj), key="i")[J(index)]
        padjDT[,i:=NULL]
        padjVals(fds, dist=i, level="site", filters=filters,
                withDimnames=FALSE) <- as.matrix(padjDT)
        
        # gene-level pval correction and FDR
        if(isTRUE(geneLevel) && 
                geneColumn %in% colnames(mcols(fds, type=type))){
            message(date(), ":  calculating gene-level pvalues ...")
            gene_pvals <- getPvalsPerGene(fds=fds, type=type, pvals=fwer_pvals,
                                            method="holm", FDRmethod=method, 
                                            geneColumn=geneColumn,
                                            BPPARAM=BPPARAM)
            pVals(fds, dist=i, level="gene", filters=filters,
                    withDimnames=FALSE) <- gene_pvals[["pvals"]]
            padjVals(fds, dist=i, level="gene", filters=filters,
                    withDimnames=FALSE) <- gene_pvals[["padj"]]
        } else if(isTRUE(geneLevel)){
            warning("Gene-level pvalues could not be calculated as column ",
                    geneColumn, "\nwas not found for the given fds object. ", 
                    "Please annotate gene symbols \nfirst using the ", 
                    "annotateRanges function.")
        }
        
        # calculate FDR for each provided subset and assign to fds
        if(!is.null(subsets)){
            stopifnot(is.list(subsets))
            stopifnot(!is.null(names(subsets)))
            for(setName in names(subsets)){
                geneListSubset <- subsets[[setName]]
                fds <- calculatePadjValuesOnSubset(fds=fds, 
                                                   genesToTest=geneListSubset,
                                                   subsetName=setName,
                                                   type=type, method=method, 
                                                   geneColumn=geneColumn, 
                                                   BPPARAM=BPPARAM)
            }
        }
    }
    
    return(fds)
}

getPvalsPerGene <- function(fds, type=currentType(fds), 
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

getSiteIndex <- function(fds, type=currentType(fds)){
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

getGeneIDs <- function(fds, type=currentType(fds), unique=TRUE, 
                        geneColumn="hgnc_symbol"){
    if(!geneColumn %in% colnames(mcols(fds, type=type))){
        stop("Did not find column '", geneColumn, "' in mcols(fds, type='", 
             type, "'). Please assign introns\nto genes first using the ",
             "annotateRanges(fds, ...) or annotateRangesWithTxDb(fds, ...) ",
             "function.")
    }
    
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

#' @describeIn FRASER This function does FDR correction only for all junctions 
#' in a certain subset of genes which can differ per sample. Requires gene 
#' symbols to have been annotated to junctions. As with the full FDR 
#' correction across all junctions, first the previously calculated 
#' junction-level p values are adjusted with Holm's method per donor or 
#' acceptor site, respectively. Then, gene-level p values are computed.
#' 
#' @param genesToTest A named list with the subset of genes to test per sample.
#'     The names must correspond to the sampleIDs in the given fds object.
#' @param subsetName The name under which the resulting FDR corrected pvalues 
#'     will be stored in metadata(fds).
#' 
#' @export
calculatePadjValuesOnSubset <- function(fds, genesToTest, subsetName, 
                                type=currentType(fds), method='BY', 
                                geneColumn="hgnc_symbol", BPPARAM=bpparam()){
    
    # check input
    stopifnot(!is.null(genesToTest))
    stopifnot(is.list(genesToTest) || is.vector(genesToTest))
    
    # replicate subset genes for each sample if given as vector
    if(!is.list(genesToTest)){
        genesToTest <- rep(list(genesToTest), ncol(fds))
        names(genesToTest) <- colnames(fds)
    }
    
    # check that names are present and correspond to samples in ods
    stopifnot(!is.null(names(genesToTest)))
    if(!all(names(genesToTest) %in% colnames(fds))){
        stop("names(genesToTest) need to be sampleIDs in the given fds object.")
    }
    
    # get genes from fds object
    fds_genes <- getGeneIDs(fds, unique=TRUE, type=type, geneColumn=geneColumn)
    ngenes <- length(fds_genes)
    
    # site index (for psi3/5)
    site_idx <- getSiteIndex(fds, type=type)
    
    # compute FDR on the given subsets of genes
    message(date(), ": starting FDR calculation on subset of genes...")
    
    # compute FDR on the given subsets of genes
    fdrSubset <- bplapply(colnames(fds), FUN=function(sampleId){
        
        # get genes to test for this sample
        genesToTestSample <- genesToTest[[sampleId]]
        padj <- rep(NA, nrow(fds))
        padj_gene <- rep(NA, ngenes)
        
        # if no genes present in the subset for this sample, return NAs
        if(is.null(genesToTestSample)){
            return(list(padj=padj, padj_gene=padj_gene))
        }
        
        # get idx of junctions corresponding to genes to test
        if(is.character(genesToTestSample)){
            rowIdx <- sort(which(fds_genes %in% genesToTestSample))
            rowIdx <- unlist(lapply(genesToTestSample, function(gene){
                idx <- which(grepl(paste0("(^|;)", gene, "(;|$)"), 
                                   mcols(fds, type=type)[, geneColumn]))
                names(idx) <- rep(gene, length(idx))
                if(length(idx) == 0 && verbose(fds) > 0){
                    warning("No introns found in fds object for gene: ", gene, 
                        " and sample: ", sampleId, ". Skipping this gene.")
                }
                return(idx)
            }))
            rowIdx <- sort(rowIdx[!duplicated(rowIdx)])
        } else{
            stop("Genes in the list to test must be a character vector ",
                 "of geneIDs.")
        }
        
        # check that rowIdx is not empty vector
        if(length(rowIdx) == 0){
            warning("No genes from the given subset found in the fds ",
                    "object for sample: ", sampleId)
            return(list(padj=padj, padj_gene=padj_gene))
        }
        
       
        
        # retrieve pvalues of introns to test
        p <- as.matrix(pVals(fds, type=type))
        if(ncol(p) == 1){
            colnames(p) <- colnames(fds)    
        }
        p <- p[rowIdx, sampleId]
        
        # FDR correction on subset
        non_dup_site_idx <- !duplicated(site_idx[rowIdx])
        padjSub <- p.adjust(p[non_dup_site_idx], method=method)
        
        # set intron FDR on subset (filled with NA for all other genes)
        padj[rowIdx] <- padjSub
        
        # gene level pvals
        dt <- data.table(sampleID=sampleId, type=type, pval=p, 
                    gene=names(rowIdx), jidx=rowIdx, site_idx=site_idx[rowIdx])
        dt <- merge(dt, 
                    data.table(site_idx=site_idx[rowIdx][non_dup_site_idx], 
                                FDR_subset=padjSub), 
                    by="site_idx")
        dt[!duplicated(dt$site_idx), 
            pval_gene:=min(p.adjust(pval, method="holm")), by="gene"]
        dt[, pval_gene := .SD[!is.na(pval_gene), unique(pval_gene)], by="gene"]
        
        # gene level FDR
        dt2 <- dt[, unique(pval_gene), by="gene"]
        dt2[, FDR_subset_gene := p.adjust(V1, method=method)]
        dt2[, gene_rowIdx := sapply(gene, function(g) which(fds_genes == g))]
        
        # set intron FDR on subset (filled with NA for all other genes)
        padj_gene[dt2[,gene_rowIdx]] <- dt2[, FDR_subset_gene]
        
        # return new FDR
        return(list(padj=padj, padj_gene=padj_gene))
        
    }, BPPARAM=BPPARAM)
    padjSub <- vapply(fdrSubset, '[[', double(nrow(fds)), 'padj')
    padjSub_gene <- vapply(fdrSubset, '[[', double(ngenes), 'padj_gene')
    
    rownames(padjSub) <- rownames(fds)
    colnames(padjSub) <- colnames(fds)
    rownames(padjSub_gene) <- fds_genes
    colnames(padjSub_gene) <- colnames(fds)
    
    # add FDR subset info to ods object and return
    padjVals(fds, type=type, level="site", subsetName=subsetName,
             withDimnames=FALSE) <- padjSub
    padjVals(fds, type=type, level="gene", subsetName=subsetName,
             withDimnames=FALSE) <- padjSub_gene
    addToAvailableFDRsubsets(fds) <- subsetName
    
    message(date(), ": finished FDR calculation on subset of genes.")
    validObject(fds)
    return(fds)
}

