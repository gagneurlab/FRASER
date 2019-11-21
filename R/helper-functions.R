#'
#' Input check functions
#'
#' Checks all user input and returns corresponding messages
#' checkFraseRDataSet
#' 
#' @return logical(1)
#' @rdname checkInputFunctions
checkFraseRDataSet <- function(fds){
    if(!is(fds, "FraseRDataSet")){
        stop("Please provide a FraseRDataSet object.")
    }
    return(invisible(TRUE))
}

#' @rdname checkInputFunctions
checkCountData <- function(fds, stop=TRUE){
    checkFraseRDataSet(fds)
    if(!all(c("rawCountsJ", "rawCountsSS") %in% assayNames(fds))){
        if(isFALSE(stop)) return(invisible(FALSE))
        stop("No counts detected! Please provide counts first.")
    }
    if(!all(paste0("rawOtherCounts_", psiTypes) %in% assayNames(fds))){
        if(isFALSE(stop)) return(invisible(FALSE))
        stop("Please compute first the total expression at each junction.")
    }
    return(invisible(TRUE))
}

#'
#' clear the files in the cache to start fresh
#'
#' @return nothing
#' @examples
#'     fds <- createTestFraseRSettings()
#'     cleanCache(fds)
#' @export
cleanCache <- function(fds, all=FALSE, cache=TRUE, assays=FALSE, results=FALSE){
    stopifnot(is(fds, "FraseRDataSet"))

    dirs2delete <- c()
    fdsDirName <- nameNoSpace(fds)
    if(cache == TRUE || all == TRUE){
        dirs2delete <- "cache"
    }
    if(assays == TRUE || all == TRUE){
        dirs2delete <- c(dirs2delete, file.path("savedObjects", fdsDirName))
    }
    if(results == TRUE || all == TRUE){
        dirs2delete <- c(dirs2delete, file.path("results", fdsDirName))
    }

    # clean cache
    for(d in dirs2delete){
        full_dir <- file.path(workingDir(fds), d)
        if(dir.exists(full_dir)){
            message(date(), ": Remove directory: '", d, "'.")
            unlink(full_dir, recursive=TRUE)
        }
    }
}

#'
#' checks if the given type is part of the correct category
#' to map correctly to a read type
#'
#' @noRd
checkReadType <- function(fds, type){

    # check if type is null or missing
    if(missing(type) | is.null(type)){
        if(verbose(fds) > 0){
            warning("Read type was not specified!",
                    "We will assume the default: 'j'")
        }
        return("j")
    }
    type <- unique(type)
    stopifnot(isScalarCharacter(type))
    correctTypes <- c(psi3="j", psi5="j", psiSite="ss")

    # check if it is already the correct type
    if(type %in% correctTypes) return(type)

    # check if psitype is given
    if(type %in% names(correctTypes)) return(correctTypes[type])

    # check assay names
    atype <- whichReadType(fds, type)
    if(!is.na(atype)) return(atype)

    # regex on the psi type
    atype <- correctTypes[vapply(names(correctTypes), FUN=grepl, type, 
                                 FUN.VALUE=logical(1))]
    if(length(atype) == 1){
        return(atype)
    }

    stop("Given read type: '", type, "' not recognized. ",
            "It needs to be 'j' (junction) or 'ss' (splice sites)",
            "\nor an existing assay name within the given object."
    )
}

#'
#' returns the corresonding PSI type to the given name
#'
#' @noRd
whichPSIType <- function(type){
    unlist(regmatches(type, gregexpr("psi(3|5|Site)", type, perl=TRUE)))
}

#'
#' returns the read type based on the given assay name
#'
#' @noRd
whichReadType <- function(fds, name){
    stopifnot(isScalarCharacter(name))

    # check writing
    if(name == "ss" | endsWith(name, "psiSite"))
        return("ss")
    if(name == "j"  | endsWith(name, "psi5") | endsWith(name, "psi3"))
        return("j")

    # check assay names
    fdsNames <- assayNames(fds)
    if(name %in% fdsNames){
        nsrNamesL <- length(assayNames(nonSplicedReads(fds)))
        fdsNamesL <- length(fdsNames)

        return(ifelse(
            which(fdsNames == name) <= fdsNamesL - nsrNamesL,
            "j",
            "ss"
        ))
    }

    stop("Could not find read type: ", name)
}

#'
#' Removes the white spaces to have a cleaner file path
#'@noRd
nameNoSpace <- function(name){
    if(is(name, "FraseRDataSet")) name <- name(name)
    stopifnot(isScalarCharacter(name))
    gsub("\\s+", "_", name, perl=TRUE)
}

#'
#' Convert to default values
#'
#' Convert all NA's of a input vector or of a
#' single dimension matrix/data.table to FALSE
#'
#' Convert NULL to NA or to another default value
#'
#' @return vector
#' @examples
#'   a <- c(TRUE, FALSE, NA, TRUE, NA)
#'   na2false(a)
#'
#'   dt <- data.table(a)
#'   na2false(dt)
#'
#'   null2na(NULL)
#'   null2na(1:10)
#'
#' @rdname na2default
#' @aliases na2false na2zero null2na null2default
#' @export
na2default <- function(x, default=FALSE){
    if(any(class(x) %in% c("DataFrame", "matrix", "data.frame"))){
        stopifnot(dim(x)[2] == 1)
        x <- as.vector(as.matrix(x)[,1])
    }
    x[is.na(x)] <- default
    return(x)
}

#' @rdname na2default
#' @export
null2default <- function(x, default=NA){
    if(is.null(x)){
        return(default)
    }
    return(x)
}

#' @rdname na2default
#' @export
na2false <- function(x){
    na2default(x, FALSE)
}

#' @rdname na2default
#' @export
na2zero <- function(x){
    na2default(x, 0)
}

#' @rdname na2default
#' @export
null2na <- function(x){
    null2default(x, NA)
}

#'
#' the qq plot function with confidence band of 5%
#' @noRd
fraserQQplotPlotly <- function(pvalues, ci=TRUE, reducePoints=FALSE,
                    sampleWise=TRUE, main="FraseR QQ-Plot"){
    if(isTRUE(reducePoints)){
        reducePoints <- c(50, 10)
    }

    # convert it to matrix if its a vector
    if(any(is(pvalues, "numeric"))){
        pvalues <- matrix(pvalues)
        colnames(pvalues) <- "observed pvalues"
    } else if(!sampleWise){
        pvalues <- t(pvalues)
    }
    if(!is.matrix(pvalues)){
        pvalues <- as.matrix(pvalues)
    }

    # convert NA to 1
    pvalues[is.na(pvalues)] <- 1

    # length of pvalues
    len_pval <- dim(pvalues)[1]

    # check colnames
    if(is.null(colnames(pvalues))){
        colnames(pvalues) <- seq_len(dim(pvalues)[2])
    }

    # my observerd and expected values
    zeroOffset <- 10e-100
    observ <- -log10(pvalues + zeroOffset)
    expect <- -log10(ppoints(len_pval))

    # create main plot object
    p <- plot_ly(type="scatter", mode="lines")

    # add theoretical trace
    p <- add_trace(p, x=expect, y=expect, mode="lines",
            line=list(color=col2hex("firebrick1")), name="theoretical-line"
    )
    p <- layout(p, title=main,
            xaxis=list(title="Expected -log<sub>10</sub>(<i>P</i>-value)"),
            yaxis=list(title="Observed -log<sub>10</sub>(<i>P</i>-value)")
    )

    # add confidence interval
    if(ci){
        if(FALSE){
            # confidence qnorm based from car::qqPlot
            o  <- abs(rnorm(10000)%%1*10^-seq(0, 3, length.out = 10000))
            o  <- sort(o)
            P  <- ppoints(o)
            n  <- length(P)
            zz <- qnorm(1-(1-0.95)/2)
            SE <- (1/dnorm(P))*sqrt(P*(1 - P)/n)
            lower <- zz*SE + P
            upper <- P*P/lower
            lower[lower==1] <- 0.999
            upper[upper==1] <- 0.999
            if(FALSE){
                plot(-log10(P), -log10(o),type="n")
                grid()
                points(-log10(P), -log10(o))
                abline(0,1, col="red")
                lines(-log10(P), -log10(upper), col="blue")
                lines(-log10(P), -log10(lower), col="green")
            }
        } else {
            # confidence qbeta based from GWASTools::qqPlot
            a <- seq_along(expect)
            upper <- -log10(qbeta(0.025, rev(a), a))
            lower <- -log10(qbeta(0.975, rev(a), a))
        }

        path  <- paste("L", c(rev(expect), expect), c(upper, rev(lower)))
        p <- layout(p, shapes=list(list(
            type="path", fillcolor="grey", opacity = 0.3,
            path=paste("M 0 0", paste(path, collapse = " "), "Z")
        )))
    }

    for(idx in seq_len(dim(pvalues)[2])){
        dat <- data.table(
            expect=expect,
            observ=sort(observ[,idx], decreasing=TRUE, na.last=TRUE)
        )
        if(length(reducePoints) > 0 & is.numeric(reducePoints)){
            nEdge <- 50
            nBy   <- 10
            if(is.numeric(reducePoints) & reducePoints[1] > 0 &
                        reducePoints[1] <= len_pval){
                nEdge <- reducePoints[1]
                if(length(reducePoints) == 2 && reducePoints[2] > 0
                            && reducePoints[2] <= 100){
                    nBy <- reducePoints[2]
                }
            }
            dat <- dat[sort(unique(c(
                seq_len(nEdge), -(nEdge-1):0+ldat, seq(1, ldat, nBy)
            )))]
        }
        p <- add_trace(p, data=dat, mode="markers",
                x=~expect, y=~observ, name=colnames(pvalues)[idx], opacity=0.3
        )
    }

    # return object
    return(p)
}

#'
#' logger function for internal use only
#' @noRd
logger <- function(type="INFO", name=flog.namespace(), ...){
    stopifnot(isScalarCharacter(type))
    type <- toupper(type)
    stopifnot(type %in% c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))

    fun <- switch(type,
            TRACE=flog.trace,
            DEBUG=flog.debug,
            INFO=flog.info,
            WARN=flog.warn,
            ERROR=flog.error,
            FATAL=flog.fatal
    )
    fun(name=name, ...)
}

#'
#' check if the given assay already exists within the object
#' @noRd
assayExists <- function(fds, assayName){
    stopifnot(isScalarCharacter(assayName))
    stopifnot(is(fds, "FraseRDataSet"))

    aexists <- assayName %in% assayNames(fds)
    if(aexists){
        message(date(), ": The ", assayName, " are already computed and will ",
                "be used now. If you want to recompute them, please remove ",
                "the following assay: ", assayName, " by issuing following ",
                "command: assays(fds)[['", assayName, "']] <- NULL")
    }
    return(aexists)
}

getAssayAsVector <- function(fds, prefix, psiType, sampleID){
    as.vector(assay(fds, paste0(prefix, psiType))[,sampleID])
}


variableJunctions <- function(fds, type, minDeltaPsi=0.1){
    psi <- K(fds, type=type)/N(fds, type=type)
    j2keep <- rowMaxs(abs(psi - rowMeans(psi, na.rm=TRUE)), na.rm=TRUE)
    j2keep >= minDeltaPsi
}

subsetKMostVariableJunctions <- function(fds, type, n){
    curX <- as.matrix(x(fds, type=type, all=TRUE, center=FALSE, 
                        noiseAlpha=NULL))
    xsd <- colSds(curX)
    nMostVarJuncs <- which(xsd >= sort(xsd, TRUE)[min(length(xsd), n*2)])
    ans <- logical(length(xsd))
    ans[sample(nMostVarJuncs, min(length(xsd), n))] <- TRUE
    ans
}

getSubsetVector <- function(fds, type, minDeltaPsi=0.1, nSubset=15000){
    # get any variable intron
    ans <- variableJunctions(fds, type, minDeltaPsi=minDeltaPsi)

    # subset most variable intron
    fds_sub <- fds[ans,,by=type]
    ans_sub <- subsetKMostVariableJunctions(fds_sub, type, nSubset)

    # set correct exclusion mask for x computation
    ans[ans] <- ans_sub
    featureExclusionMask(fds) <- exMask
}

pasteTable <- function(x, ...){
    tab <- table(x, ...)
    paste(names(tab), tab, collapse="\t", sep=": ")
}

#'
#' Map between individual seq level style and dataset common one
#' for counting and aggregating the reads
#' @noRd
checkSeqLevelStyle <- function(gr, fds, sampleID, sampleSpecific=FALSE){
    if(!"SeqLevelStyle" %in% colnames(colData(fds))){
        return(gr)
    }
    style <- colData(fds)[sampleID,"SeqLevelStyle"]
    if(isFALSE(sampleSpecific)){
        style <- names(sort(table(colData(fds)[,"SeqLevelStyle"]), TRUE)[1])
        if(length(unique(colData(fds)[,"SeqLevelStyle"])) > 1){
            gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
        }
    }

    seqlevelsStyle(gr) <- style
    gr
}

uniformSeqInfo <- function(grls){
    tmpSeqlevels <- unique(data.table(
        seqlevel  = unlist(lapply(grls, seqlevels)),
        seqlength = unlist(lapply(grls, seqlengths))
    )[order(seqlevel)])

    if(any(duplicated(tmpSeqlevels[,seqlevel]))){
        stop("There are non uniq chromosomes in this dataset!")
    }

    ans <- lapply(grls, function(x){
        seqlevels(x)  <- tmpSeqlevels[,seqlevel]
        seqlengths(x) <- tmpSeqlevels[,seqlength]
        x
    })
    ans
}

getHDF5ChunkSize <- function(fds, assayName){
    h5obj <- H5Fopen(getFraseRHDF5File(fds, assayName), flags="H5F_ACC_RDONLY")
    ans <- rhdf5:::H5Dchunk_dims(h5obj&assayName)
    H5Fclose(h5obj)
    ans
}

getMaxChunks2Read <- function(fds, assayName, max=15, axis=c("col", "row")){
    if(!any(c("DelayedArray", "DelayedMatrix") %in%
            class(assay(fds, assayName)))){
        if(axis == "col"){
            return(ceiling(ncol(assay(fds, assayName))/bpnworkers(bpparam())))
        }
        return(ceiling(nrow(assay(fds, assayName))/bpnworkers(bpparam())))
    }

    axis <- match.arg(axis)
    dims <- getHDF5ChunkSize(fds, assayName)
    if(axis == "col"){
        ans <- dims[2]
    } else {
        ans <- dims[1]
    }
    max(1, ans/ceiling(ans/max))
}

getSamplesByChunk <- function(fds, sampleIDs, chunkSize){
    chunks <- trunc(0:(ncol(fds)-1)/chunkSize)
    ans <- lapply(0:max(chunks), function(x){
        intersect(sampleIDs, samples(fds)[chunks == x])
    })
    ans[vapply(ans, length, integer(1)) >0]
}

checkNaAndRange <- function(x, min=-Inf, max=Inf, scalar=TRUE, na.ok=FALSE){
    xname <- deparse(substitute(x))
    if(isTRUE(scalar) & !isScalarValue(x)){
        stop(xname, " should be a scalar value!")
    }
    if(any(is.na(x)) && isFALSE(na.ok)){
        stop(xname, " contains NA values, which is not allowed.")
    }
    if(sum(!is.na(x)) == 0){
        return(invisible(TRUE))
    }
    if(!is.numeric(x[!is.na(x)])){
        stop(xname, " should be numeric!")
    }
    if(x[!is.na(x)] < min){
        stop(xname, " should be bigger than ", min)
    }
    if(x[!is.na(x)] > max){
        stop(xname, " should be smaller than ", max)
    }
    invisible(TRUE)
}

putCounts2Memory <- function(fds, type=currentType(fds)){
    counts(fds, type=type, side="other", HDF5=FALSE) <-
            as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- 
            as.matrix(counts(fds, type=type, side="ofInterest"))
    fds
}


plotBasePlot <- function(ggplot, basePlot=FALSE){
    if(isFALSE(basePlot)){
        if(!require(plotly)){
            stop("Please install plotly, if you use the option basePlot=FALSE!")
        }
        ggplot$labels <- lapply(ggplot$labels, function(x){
                if(typeof(x) == "expression"){
                    warning("Found expression for plotly. Please adapt it!")
                    return(as.character(x))
                }
                x})
        return(plotly::ggplotly(ggplot, tooltip="text"))
    }
    ggplot
}
