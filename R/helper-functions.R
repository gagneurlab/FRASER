########
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## Helperfunctions to convert or extract data within the FraseR package
##

#'
#' clear the files in the cache to start fresh
#'
#' @examples
#'     fds <- createTestFraseRSettings()
#'     cleanCache(fds)
#'
#' @export
cleanCache <- function(fds, all=FALSE, cache=TRUE, assays=FALSE, results=FALSE){
    stopifnot(class(fds) == "FraseRDataSet")

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
    if(!is.na(atype)) return(atype)

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
#' convert all NA's of a input vector or of a
#' single dimension matrix/data.table to FALSE
#'
#' @examples
#'   a <- c(TRUE, FALSE, NA, TRUE, NA)
#'   na2false(a)
#'
#'   dt <- data.table(a)
#'   na2false(dt)
#'
#' @export
na2false <- function(x){
    na2default(x, FALSE)
}
na2zero <- function(x){
    na2default(x, 0)
}
na2default <- function(x, default=FALSE){
    if(any(class(x) %in% c("DataFrame", "matrix", "data.frame"))){
        stopifnot(dim(x)[2] == 1)
        x <- as.vector(as.matrix(x)[,1])
    }
    x[is.na(x)] <- default
    return(x)
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
    if(any(class(pvalues) == "numeric")){
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
        colnames(pvalues) <- 1:dim(pvalues)[2]
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
            a <- 1:length(expect)
            upper <- -log10(qbeta(0.025, rev(a), a))
            lower <- -log10(qbeta(0.975, rev(a), a))
        }

        path  <- paste("L", c(rev(expect), expect), c(upper, rev(lower)))
        p <- layout(p, shapes=list(list(
            type="path", fillcolor="grey", opacity = 0.3,
            path=paste("M 0 0", paste(path, collapse = " "), "Z")
        )))
    }

    for(idx in 1:dim(pvalues)[2]){
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
                1:nEdge, -(nEdge-1):0+ldat, seq(1, ldat, nBy)
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
#'
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
#'
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
