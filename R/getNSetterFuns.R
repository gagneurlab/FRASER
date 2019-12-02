#' Getter/Setter functions
#' 
#' This is a collection of small accessor/setter functions for easy access to
#' the values within the FraseR model.
#' 
#' @param fds An FraseRDataSet object.
#' @param type The type of psi (psi5, psi3 or psiSite)
#' @param byGroup If TRUE, aggregation by donor/acceptor site will be done.
#' @param dist Distribution for which the p-values should be extracted.
#' @param value The new value to be assigned. 
#' @return A (delayed) matrix or vector dependent on the type of data retrieved.
#' 
#' @name getter_setter_functions
#' @rdname getter_setter_functions
#' @aliases zScore, pVals, padjVals, rho, bestQ
#' 
#' @examples 
#' fds <- createTestFraseRDataSet()
#' dontWriteHDF5(fds)
#' dontWriteHDF5 <- TRUE
#' currentType(fds) <- "psi5"
#' 
#' bestQ(fds)
#' rho(fds)
#' # get statistics
#' pVals(fds)
#' padjVals(fds)
#' zScores(fds)
#' 
NULL

#'
#' Feature exclusion
#'
#' To remove certain junctions from being used in the train step of the
#' encoding dimension we can set the \code{featureExclusion} vector to
#' \code{FALSE}. This can be helpfull if we have local linkage between
#' features which we do not want to model by the autoencoder.
#'
#' @param fds A FraseRDataSet object
#' @param type The psi type.
#' @param value A logical vector of the length of the features. If
#'             \code{TRUE}, the corresponding feature will be excluded
#'             from the encoding dimension fit.
#' @return The exclusion vector
#'
#' @rdname featureExclusionMask
#' @aliases featureExclusionMask, `featureExclusionMask<-`
#'
#' @examples
#' fds <- makeExampleFraseRDataSet()
#' featureExclusionMask(fds) <- sample(c(FALSE, TRUE), nrow(fds), replace=TRUE)
#'
#' featureExclusionMask(fds)
#'
#' @export featureExclusionMask
featureExclusionMask <- function(fds, type=currentType(fds)){
    ans <- rep(TRUE, nrow(mcols(fds, type=type)))
    if(paste0('featureExclude_', type) %in% colnames(mcols(fds, type=type))){
        ans <- mcols(fds, type=type)[[paste0('featureExclude_', type)]]
    }
    # TODO names(ans, type=type) <- rownames(fds, type=type)
    return(ans)
}

#' @rdname featureExclusionMask
#' @export "featureExclusionMask<-"
`featureExclusionMask<-` <- function(fds, type=currentType(fds), value){
    if(isScalarLogical(value)){
        value <- rep(value, nrow(mcols(fds, type=type)))
    }
    mcols(fds, type=type)[[paste0('featureExclude_', type)]] <- value
    return(fds)
}

#' @rdname counts
#' @export
K <- function(fds, type=currentType(fds)){
    K <- counts(fds, type=type, side="ofInterest")
    return(K);
}

#' @rdname counts
#' @export
N <- function(fds, type=currentType(fds)){
    N <- K(fds, type=type) + counts(fds, type=type, side="other")
    return(N);
}

x <- function(fds, type=currentType(fds), all=FALSE,
                    noiseAlpha=currentNoiseAlpha(fds), center=TRUE){
    K <- K(fds, type=type)
    N <- N(fds, type=type)

    # compute logit ratio with pseudocounts
    x <- as.matrix(t((K + pseudocount())/(N + (2*pseudocount()))))
    x <- qlogis(x)
    
    if(any(is.infinite(x))){
        x[is.infinite(x) & x > 0] <- NA
        x[is.na(x)] <- max(x, na.rm=TRUE) + 1
    }

    # corrupt x if required
    if(!is.null(noiseAlpha)){
        noise <- noise(fds, type=type)
        if(is.null(noise)){
            noise <- matrix(rnorm(prod(dim(x))), ncol=ncol(x), nrow=nrow(x))
            noise(fds, type=type) <- noise
        }
        x <- x + t(colSds(x) * noiseAlpha * t(noise))
    }

    if(isFALSE(all)){
        x <- x[,featureExclusionMask(fds, type=type)]
    }
    if(isTRUE(center)){
        x <- t(t(x) - colMeans2(x))
    }

    return(x)
}

H <- function(fds, type=currentType(fds), noiseAlpha=NULL){
    x(fds, all=FALSE, type=type, noiseAlpha=noiseAlpha) %*% E(fds, type=type)
}

`D<-` <- function(fds, type=currentType(fds), value){
    if(!is.matrix(value)){
        value <- matrix(value, nrow=nrow(fds))
    }
    metadata(fds)[[paste0('D_', type)]] <- value
    return(fds)
}

D <- function(fds, type=currentType(fds)){
    return(metadata(fds)[[paste0('D_', type)]])
}

`E<-` <- function(fds, type=currentType(fds), value){
    if(!is.matrix(value)){
        value <- matrix(value, nrow=sum(featureExclusionMask(fds, type=type)))
    }
    metadata(fds)[[paste0('E_', type)]] <- value
    return(fds)
}

E <- function(fds, type=currentType(fds)){
    return(metadata(fds)[[paste0('E_', type)]])
}

`b<-` <- function(fds, type=currentType(fds), value){
    mcols(fds, type=type)[[paste0('b_', type)]] <- value
    return(fds)
}

b <- function(fds, type=currentType(fds)){
    return(mcols(fds, type=type)[[paste0('b_', type)]])
}

`rho<-` <- function(fds, type=currentType(fds), value){
    mcols(fds, type=type)[[paste0('rho_', type)]] <- value
    return(fds)
}

#' @describeIn getter_setter_functions Returns the fitted rho values for the 
#' beta-binomial distribution
#' @export
rho <- function(fds, type=currentType(fds)){
    return(mcols(fds, type=type)[[paste0('rho_', type)]])
}

predictMu <- function(fds, type=currentType(fds), noiseAlpha=NULL){
    y <- predictY(fds, type=type, noiseAlpha=noiseAlpha)
    mu <- predictMuCpp(y)
    return(t(mu))
}

predictY <- function(fds, type=currentType(fds), noiseAlpha=NULL){
    D <- D(fds, type=type)
    b <- b(fds, type=type)
    H <- H(fds, type=type, noiseAlpha=noiseAlpha)

    y <- predictYCpp(as.matrix(H), D, b)

    return(t(y))
}


`setAssayMatrix<-` <- function(fds, name, type, ..., value){
    if(!is.matrix(value)){
        value <- matrix(value, ncol=ncol(fds), nrow=nrow(mcols(fds, type=type)))
    }
    if(is.null(colnames(value))){
        colnames(value) <- colnames(fds)
    }
    if(is.null(rownames(value))){
        rownames(value) <- rownames(counts(fds, type=type))
    }
    if(missing(name)){
        name <- type
    } else {
        name <- paste(name, type, sep="_")
    }
    assay(fds, name, ...) <- value
    fds
}

getAssayMatrix <- function(fds, name, type, byGroup=FALSE){
    if(missing(name)){
        name <- type
    } else {
        name <- paste(name, type, sep="_")
    }
    ans <- assay(fds, name)
    idx <- TRUE
    if(isTRUE(byGroup)){
        idx <- !duplicated(getSiteIndex(fds, type=type))
    }
    ans[idx,]
}

#' @describeIn getter_setter_functions This returns the calculated z-scores.
#' @export
zScores <- function(fds, type=currentType(fds), byGroup=FALSE){
    ans <- getAssayMatrix(fds, name='zScores', type=type)
    if(isTRUE(byGroup)){
        ans <- getAbsMaxByGroup(fds, type, ans)
    }
    ans
}

`zScores<-` <- function(fds, type=currentType(fds), ..., value){
    setAssayMatrix(fds, name="zScores", type=type, ...) <- value
    return(fds)
}
#' @describeIn getter_setter_functions This returns the calculated p-values.
#' @export
pVals <- function(fds, type=currentType(fds), 
                    dist="BetaBinomial", byGroup=FALSE){
    dist <- match.arg(dist, choices=c("BetaBinomial", "Binomial", "Normal"))
    getAssayMatrix(fds, paste0("pvalues", dist), type=type, byGroup=byGroup)
}

`pVals<-` <- function(fds, type=currentType(fds),
                    dist="BetaBinomial", ..., value){
    dist <- match.arg(dist, choices=c("BetaBinomial", "Binomial", "Normal"))
    setAssayMatrix(fds, name=paste0("pvalues", dist), type=type, ...) <- value
    return(fds)
}

#' @describeIn getter_setter_functions This returns the adjusted p-values.
#' @export
padjVals <- function(fds, type=currentType(fds),
                    dist=c("BetaBinomial"), byGroup=FALSE){
    dist <- match.arg(dist, choices=c("BetaBinomial", "Binomial", "Normal"))
    getAssayMatrix(fds, paste0("pajd", dist), type=type, byGroup=byGroup)
}

`padjVals<-` <- function(fds, type=currentType(fds),
                    dist="BetaBinomial", ..., value){
    dist <- match.arg(dist, choices=c("BetaBinomial", "Binomial", "Normal"))
    setAssayMatrix(fds, name=paste0("pajd", dist), type=type, ...) <- value
    return(fds)
}

#' @describeIn getter_setter_functions This returns the fitted mu (i.e. psi) 
#' values.
#' @export
predictedMeans <- function(fds, type=currentType(fds)){
    return(getAssayMatrix(fds, name="predictedMeans", type=type))
}

`predictedMeans<-` <- function(fds, type=currentType(fds), ..., value){
    setAssayMatrix(fds, name="predictedMeans", type=type, ...) <- value
    return(fds)
}

#' @describeIn getter_setter_functions Returns the difference between the 
#' observed and the fitted psi values.
#' @export
deltaPsiValue <- function(fds, type=currentType(fds)){
    return(assay(fds, type) - predictedMeans(fds, type=type))
}


#'
#' Set/get psi type
#'
#' Set and returns the psi type that is to be used within several methods in 
#' the FraseR package.
#'
#' @param fds FraseRDataSet
#' @param value If given, is sets the psi type. It is required to be a scalar 
#' character.
#' @return character() (getter) or FraseRDataSet(setter)
#' @examples
#'   # get data
#'   fds <- makeSimulatedFraserDataSet(m=50, j=1000)
#' 
#'   # set
#'   currentType(fds) <- "psi5"
#'
#'   # get
#'   currentType(fds)
#'
#' @export
currentType <- function(fds){
    return(metadata(fds)[['currentType']])
}

#' @rdname currentType
#' @export 
`currentType<-` <- function(fds, value){
    stopifnot(isScalarCharacter(whichPSIType(value)))
    metadata(fds)[['currentType']] <- whichPSIType(value)
    return(fds)
}

#'
#' Set/get global pseudo count option
#'
#' Set and returns the pseudo count used within the FraseR fitting procedure.
#'
#' @param value If given, is sets the psuedocount. It requires a positive
#'                number and rounds it to the next integer.
#' @return numeric scalar
#' @examples
#' # set
#' pseudocount(4L)
#'
#' # get
#' pseudocount()
#'
#' @export
pseudocount <- function(value=NULL){
    # return if not provided
    if(is.null(value)){
        ans <- options()[['FraseR.pseudoCount']]
        if(isScalarNumeric(ans)){
            return(ans)
        }
        return(1)
    }

    # set pseudo count if provided
    stopifnot(isScalarNumeric(value))
    stopifnot(value >= 0)
    value <- as.integer(value)
    options('FraseR.pseudoCount'=value)
    devNULL <- .setPseudoCount(value)
    stopifnot(value == devNULL)

    invisible(value)
}

currentNoiseAlpha <- function(fds){
    return(metadata(fds)[['noiseAlpha']])
}

`currentNoiseAlpha<-` <- function(fds, value){
    metadata(fds)[['noiseAlpha']] <- value
    return(fds)
}

noise <- function(fds, type=currentType(fds)){
    return(t(getAssayMatrix(fds, name="noise", type=type)))
}

`noise<-` <- function(fds, type=currentType(fds), HDF5=FALSE, ..., value){
    if(!is.matrix(value)){
        value <- matrix(value, nrow=nrow(mcols(fds, type=type)), ncol=ncol(fds))
    }
    setAssayMatrix(fds, name='noise', type=type, HDF5=HDF5) <- t(value)
    return(fds)
}


hyperParams <- function(fds, type=currentType(fds), all=FALSE){
    ans <- metadata(fds)[[paste0("hyperParams_", type)]]
    if(is.null(ans)){
        return(ans)
    }
    if(isFALSE(all)){
        ans <- ans[aroc == max(aroc)][1]
    }
    ans
}

`hyperParams<-` <- function(fds, type=currentType(fds), value){
    metadata(fds)[[paste0("hyperParams_", type)]] <- value
    return(fds)
}

#' @describeIn getter_setter_functions This returns the optimal size of the 
#' latent space according to the hyperparameter optimization or a simple 
#' estimate of about a tenth of the number of samples if the hyperparameter 
#' opimization was not run yet.
#' @export
bestQ <- function(fds, type=currentType(fds)){
    ans <- hyperParams(fds, type=type)[1,q]
    if(is.null(ans) || is.na(ans)){
        warnings("Please set q by estimating it correctly.")
        ans <- min(100, max(2, round(ncol(fds)/10)))
    }
    return(as.integer(ans))
}

bestNoise <- function(fds, type=currentType(fds)){
    ans <- hyperParams(fds, type=type)[1,noise]
    if(is.null(ans)){
        warnings("Please set noise by estimating it correctly.")
        ans <- 1
    }
    as.numeric(as.character(ans))
}

#' @describeIn getter_setter_functions Gets the current value of whether the 
#' assays should be stored as hdf5 files.
#' @export
dontWriteHDF5 <- function(fds){
    return(metadata(fds)[['dontWriteHDF5']])
}

#' @describeIn getter_setter_functions Sets whether the assays should be stored 
#' as hdf5 files.
#' @export
`dontWriteHDF5<-` <- function(fds, value){
    metadata(fds)[['dontWriteHDF5']] <- isTRUE(value)
    return(fds)
}

getTrueOutliers <- function(fds, type, BPPARAM=bpparam(), byGroup=FALSE){
    ans <- getAssayMatrix(fds, "trueOutliers", type)
    if(isTRUE(byGroup)){
        ans <- getAbsMaxByGroup(fds, type, ans, BPPARAM)
    }
    
    # remove secondary injections -> -1/0/+1 instead of -2/-1/0/+1/+2
    pmin(pmax(ans, -1), 1)
}

getTrueDeltaPsi <- function(fds, type, BPPARAM=bpparam(), byGroup=FALSE){
    ans <- getAssayMatrix(fds, "trueDeltaPSI", type)
    if(isTRUE(byGroup)){
        ans <- getAbsMaxByGroup(fds, type, ans, BPPARAM)
    }
    ans
}

getAbsMaxByGroup <- function(fds, type, mat, index=NULL, BPPARAM=bpparam()){
    if(is.null(index)){
        index <- getSiteIndex(fds, type)
    }

    dt <- cbind(data.table(id=index), as.data.table(mat))
    values <- matrix(ncol=ncol(mat), unlist(bplapply(colnames(mat), 
            BPPARAM=BPPARAM,
            function(i){
                    dttmp <- dt[,.(.I, id, value=get(i), abs=abs(get(i)))]
                    dttmp[,maxVal:=value[abs == max(abs)][1], by=id]
                    dttmp[!duplicated(id)][order(I)][,value]
            })))
    
    colnames(values) <- colnames(mat)
    rownames(values) <- index[!duplicated(index)]
    return(values)
}

getByGroup <- function(fds, type, value){
    index <- getSiteIndex(fds, type)
    idx   <- !duplicated(index)
    return(value[idx,])
}

getDeltaPsi <- function(fds, type, byGroup=FALSE){
    mu <- predictedMeans(fds, type)
    dataPsi <- (K(fds, type) + pseudocount())/(N(fds, type) + 2*pseudocount())
    deltaPSI <- dataPsi - mu
    if(isTRUE(byGroup)){
        deltaPSI <- getAbsMaxByGroup(fds, psiType, deltaPSI)
    }
    deltaPSI
}


# calculate FraseR weights
calcFraseRWeights <- function(fds, psiType){
    k <- as.matrix(K(fds, psiType))
    n <- as.matrix(N(fds, psiType))
    mu <- t(predictMu(fds, psiType))
    rho <- rho(fds, psiType)
    dataPsi <- plogis(t(
            x(fds, type=psiType, all=TRUE, center=FALSE, noiseAlpha=NULL)))

    # pearson residuals for BB
    # on counts of success k
    # r <- ((k+pseudocount()) - (n+2*pseudocount()) * mu) / sqrt(
    #       (n+2*pseudocount()) * mu * (1-mu) *
    #       (1+((n+2*pseudocount())-1)*rho))
    # on probability of success mu
    r <- (dataPsi - mu) / sqrt(
            mu * (1-mu) * (1+((n+2*pseudocount())-1)*rho) /
            (n+2*pseudocount()))

    # weights according to Huber function (as in edgeR)
    c <- 1.345; # constant, as suggested in edgeR paper
    w <- ifelse(abs(r) > c, c/abs(r) , 1)

    return(w)
}

# get FraseR weights
weights <- function(fds, type){
    return(getAssayMatrix(fds, "weights", type))
}

# set FraseR weights
`weights<-` <- function(fds, type=currentType(fds), ..., value){
    setAssayMatrix(fds, name="weights", type=type, ...) <- value
    return(fds)
}

getIndexFromResultTable <- function(fds, resultTable, padj.method="holm"){
    type <- as.character(resultTable$type)
    target <- makeGRangesFromDataFrame(resultTable)
    if(type == "psiSite"){
        gr <- granges(asSE(nonSplicedReads(fds)))
    } else {
        gr <- granges(asSE(fds))
    }

    hits <- findOverlaps(target, gr, type="equal")
    ov <- to(hits)
    if(!isScalarInteger(ov)){
        stop("Can not find the given range within the FraseR object.")
    }
    ov
}

getPlottingDT <- function(fds, axis=c("row", "col"), type=NULL, result=NULL,
                    idx=NULL, aggregate=FALSE, BPPARAM=SerialParam(),
                    Ncpus=min(10, getDTthreads()), ...){
    if(!is.null(result)){
        type <- as.character(result$type)
        idx  <- getIndexFromResultTable(fds, result)
    }

    axis <- match.arg(axis)
    idxrow <- idx
    idxcol <- TRUE
    if(axis == "col"){
        idxcol <- idx
        if(is.character(idx)){
            idxcol <- colnames(fds) %in% idx
        }
        idxrow <- TRUE
    }

    k <- K(fds, type)[idxrow, idxcol]
    n <- N(fds, type)[idxrow, idxcol]
    
    spliceID <- getSiteIndex(fds, type=type)[idxrow]
    feature_names <- rownames(mcols(fds, type=type))[idxrow]
    if("hgnc_symbol" %in% colnames(mcols(fds, type=type))){
        feature_names <- mcols(fds, type=type)[idxrow,"hgnc_symbol"]
    }
    if(is.null(feature_names)){
        feature_names <- as.character(seq_row(mcols(fds, type=type)))[idxrow]
    }

    dt <- data.table(
            idx       = idx,
            sampleID  = factor(as.character(colnames(K(fds, type))[idxcol])),
            spliceID  = factor(spliceID),
            featureID = factor(feature_names),
            type      = factor(type),
            k         = c(k),
            n         = c(n),
            pval      = c(pVals(fds, type=type)[idxrow, idxcol]),
            padj      = c(padjVals(fds, type=type)[idxrow, idxcol]),
            zscore    = c(zScores(fds, type=type)[idxrow, idxcol]),
            obsPsi    = c((k + pseudocount())/(n + 2*pseudocount())),
            predPsi   = c(predictedMeans(fds, type)[idxrow, idxcol]))
    dt[, deltaPsi:=obsPsi - predPsi]

    # if requested return gene p values (correct for multiple testing again)
    if(isTRUE(aggregate)){
        dt <- dt[!is.na(featureID)]

        # correct by gene and take the smallest p value
        dt <- rbindlist(mclapply(unique(dt[,sampleID]), mc.cores=Ncpus,
            FUN=function(x){
                    dttmp <- dt[sampleID == x]
                    dttmp[, pval:=p.adjust(pval, method="holm"),
                            by="sampleID,featureID,type"]
                    dttmp <- dttmp[order(sampleID, featureID, type, pval)][
                            !duplicated(data.table(sampleID, featureID, type))]
                    dttmp <- dttmp[, padj:=p.adjust(pval, method="BY"),
                            by="sampleID,type"]
                    dttmp
            }))
    }

    # add aberrant information to it
    aberrantVec <- aberrant(fds, ..., padjVals=dt[,.(padj)],
        dPsi=dt[,.(deltaPsi)], zscores=dt[,.(zscore)])
    dt[,aberrant:=aberrantVec]

    # return object
    dt
}


#'
#' Verbosity level of the FraseR package
#'
#' Dependend on the level of verbosity the algorithm reports more or less to
#' the user. 0 means being quiet and 10 means everything.
#'
#' @param fds FraseRDataSet
#' @param value The level of verbosity, between 0 and 10. TRUE/FALSE are also 
#' accepted.
#'
#' @rdname verbose
#' @return numeric(1) (getter) or FraseRDataSet (setter)
#' @examples
#' fds <- createTestFraseRSettings()
#' verbose(fds)
#'
#' verbose(fds) <- 2
#' verbose(fds)
#'
#' @export
verbose <- function(fds){
    if("verbosity" %in% names(metadata(fds))){
        return(metadata(fds)[["verbosity"]])
    }
    return(0)
}

#' @rdname verbose
#' @export
`verbose<-` <- function(fds, value){
    verbose <- value
    if(is.logical(verbose)){
        verbose <- verbose + 0
    }
    checkNaAndRange(verbose, min=0, max=10, na.ok=FALSE)
    metadata(fds)[["verbosity"]] <- floor(verbose)
    return(fds)
}
