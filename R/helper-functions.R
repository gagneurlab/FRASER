########
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## Helperfunctions to convert or extract data within the FraseR package
##

#'
#' clear the files in the cache to start fresh
#' @export
cleanCache <- function(fds, ...){
    stopifnot(class(fds) == "FraseRDataSet")
    # clean cache
    cacheDir <- file.path(workingDir(fds), "cache")
    if(dir.exists(cacheDir)){
        unlink(cacheDir, recursive=TRUE)
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
        warning("Read type was not specified! We will assume the default: 'j'")
        return("j")
    }

    stopifnot(isScalarCharacter(type))
    correctTypes <- c(psi3="j", psi5="j", psiSite="ss")

    # check if it is already the correct type
    if(type %in% correctTypes) return(type)

    # check if psitype is given
    if(type %in% names(correctTypes)) return(correctTypes[type])

    # check assay names
    atype <- whichReadType(fds, type)
    if(isCorrectType(atype)) return(atype)

    stop("Given read type: '", type, "' not recognized. ",
            "It needs to be 'j' (junction) or 'ss' (splice sites)",
            "\nor an existing assay name within the given object."
    )
}

#'
#' returns the read type based on the given assay name
#'
#' @noRd
whichReadType <- function(fds, name){
    stopifnot(isScalarCharacter(name))
    fdsNames <- assayNames(fds)
    if(!name %in% fdsNames){
        return(NA)
    }
    nsrNamesL <- length(assayNames(nonSplicedReads(fds)))
    fdsNamesL <- length(fdsNames)

    return(ifelse(
        which(fdsNames == name) <= fdsNamesL - nsrNamesL,
        "j",
        "ss"
    ))
}

#'
#' Removes the white spaces to have a cleaner file path
#'@noRd
nameNoSpace <- function(name){
    if(class(name) == "FraseRDataSet") name <- name(name)
    stopifnot(isScalarCharacter(name))
    gsub("\\s+", "_", name, perl=TRUE)
}

#'
#' convert the input of NA to FALSE
#'
#' @export
na2false <- function(x){
    if(any(class(x) %in% c("DataFrame", "matrix", "data.frame"))){
        stopifnot(dim(x)[2] == 1)
        x <- as.vector(as.matrix(x)[,1])
    }
    x[is.na(x)] <- FALSE
    return(x)
}

#'
#' the qq plot function
#' @noRd
fraserQQplotPlotly <- function(pvalues, zscores=NULL, zscoreCutoff=0,
                reducePoints=0, main="FraseR QQ-Plot", ...){
    # cutoff by zscore before generating the qq plot
    lpval <- ifelse(any(class(pvalues) == "numeric"),
            length(pvalues), dim(pvalues)[1]
    )
    goodZscores <- rep(TRUE, lpval)
    if(!is.null(zscores)){
        stopifnot(zscoreCutoff >= 0 && zscoreCutoff <= 100)
        absZscores <- as.matrix(abs(zscores))
        goodZscores <- na2false(rowMaxs(absZscores, na.rm=TRUE) > zscoreCutoff)
    }

    # my observerd and expected values
    zeroOffset <- 10e-100
    observ <- -log10(pvalues[goodZscores] + zeroOffset)
    expect <- -log10(ppoints(sum(goodZscores)))

    p <- plot_ly(type="scattergl", mode="lines")
    for(s in colnames(pvalues)){
        dat <- data.table(
            expect=expect,
            observ=sort(observ[,get(s)], decreasing=TRUE, na.last=TRUE)
        )
        dat <- dat[!is.na(observ)]
        ldat <- nrow(dat)
        if(reducePoints){
            nEdge <- 50
            nBy   <- 10
            if(numeric(reducePoints) && reducePoints[1] > 0
                        && reducePoints[1] <= 100){
                nEdge <- reducePoints[1]
                if(length(reducePoints) == 2 && reducePoints[2] > 0
                            && reducePoints[2] <= 100){
                    nBy <- reducePoints[2]
                }
            }
            dat <- dat[unique(c(1:nEdge, -(nEdge-1):0+ldat, seq(1, ldat, nBy)))]
        }
        p <- add_trace(p, data=dat, mode="markers",
                       x=~expect, y=~observ, name=s, opacity=0.3
        )
    }
    p <- add_trace(p, x=expect, y=expect, mode="lines", name="theoretical-line")
    p <- layout(p, title=main,
        xaxis=list(title="Expected -log10(<i>P-value)"),
        yaxis=list(title="Observed -log10(<i>P-value)")
    )
    p
    return(p)
}


