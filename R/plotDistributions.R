#'
#' Visualisation of the distributions of the statistics used within FraseR.
#'
#' This can be used for quality control plots as well.
#'
#' @param fds A FraseR object
#' @param gr  The range of interest to visualize.
#'              Needs to match a junction or splice site within the FraseR obj.
#'              Also it needs to have a type column specifying
#'              the site of interest ("psi3", "psi5", "psiSite").
#' @param sample A vector of sampleIDs which should be highlighted.
#'              If this value is not set but the gr object contains a sampleID
#'              column it will take this column to highlight samples.
#' @param ... Additional paramerters which are pased on to
#'              the underlying plotting functions
#'
#' @export
plotJunctionDistribution <- function(fds, gr, type=gr$type, sampleIDs=NULL,
            rmZeroCts=FALSE, valueVsCounts=FALSE, qqplot=FALSE){
    stopifnot(is(fds, "FraseRDataSet"))
    if(is.data.table(gr)){
        gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE)
    }
    stopifnot(is(gr, "GRanges"))
    stopifnot(length(gr) == 1)

    if(is.null(sampleIDs) & "sampleID" %in% colnames(mcols(gr))){
        sampleIDs <-  gr$sampleID
    }

    old.par <- par(no.readonly=TRUE)

    data <- getPlotDistributionData(gr, fds, type, rmZeroCts)

    par(mfrow=c(2,ifelse(qqplot, 3, 2)), cex=1)
    plotSampleRank(gr, fds, type=type, sample=sampleIDs, rmZeroCts=rmZeroCts,
            data=data)
    plotSampleRank(gr, fds, type=type, sample=sampleIDs, rmZeroCts=rmZeroCts,
            delta=TRUE, data=data)
    plotCountsAtSite(gr, fds, type=type, sample=sampleIDs, data=data)
    if(!missing(valueVsCounts) & (any(valueVsCounts %in% assayNames(fds)) |
                    is.logical(valueVsCounts) & valueVsCounts != FALSE)){
        type2plot <- ifelse(isTRUE(valueVsCounts), type, valueVsCounts)
        plotValueVsCounts(gr, fds, type2plot, sampleIDs, plotlog=TRUE,
                rmZeroCts=rmZeroCts, data=data)
    } else {
        plotSampleRank(gr, fds, type=paste0("zscore_", type), sample=sampleIDs,
                rmZeroCts=rmZeroCts, data=data)
    }
    if(qqplot){
        plotQQplot(gr, fds, type=type, data=data)
    }

    invisible(eval(par(old.par), envir = parent.frame()))
}

#'
#' get data for plot
#' @noRd
getPlotDistributionData <- function(gr, fds, type, rmZeroCts=FALSE){
    se <- as(fds, "RangedSummarizedExperiment")
    if(type %in% "psiSite") {
        se <- nonSplicedReads(fds)
    }

    rctsname     <- paste0("rawCounts", toupper(checkReadType(fds, type)))
    roctsname    <- paste0("rawOtherCounts_", whichPSIType(type))
    pvaluename   <- paste0("pvalue_", type)
    psiname      <- type
    deltapsiname <- paste0("delta_", type)
    zscorename   <- paste0("zscore_", type)

    mapping      <- c("pvalues", "psi", "deltaPsi", "zscore", "rcts", "rocts")
    names(mapping) <- c(pvaluename, psiname, deltapsiname, zscorename, rctsname,
            roctsname)

    ov <- from(findOverlaps(granges(se), gr, type="equal"))
    if(!isScalarInteger(ov)){
        stop("Can not find the given GRange object within the FraseR object.")
    }
    seOfInterest <- se[ov]

    rcts     = as.matrix(assays(seOfInterest)[[rctsname]])[1,]
    rocts    = as.matrix(assays(seOfInterest)[[roctsname]])[1,]
    if(rmZeroCts){
        se <- se[,rcts + rocts > 0]
    }

    return(list(
        se       = seOfInterest,
        mapping  = mapping,
        pvalues  = as.matrix(assays(seOfInterest)[[pvaluename]])[1,],
        psi      = as.matrix(assays(seOfInterest)[[psiname]])[1,],
        deltaPsi = as.matrix(assays(seOfInterest)[[deltapsiname]])[1,],
        zscore   = as.matrix(assays(seOfInterest)[[zscorename]])[1,],
        rcts     = as.matrix(assays(seOfInterest)[[rctsname]])[1,],
        rocts    = as.matrix(assays(seOfInterest)[[roctsname]])[1,],
        alpha    = mcols(seOfInterest)[,paste0(type, "_alpha")],
        beta     = mcols(seOfInterest)[,paste0(type, "_beta")]
    ))
}


#'
#' plot count distribution
#'
#' @noRd
plotCountsAtSite <- function(gr, fds, type, sample=NULL, plotLog=TRUE,
                    plotLegend=TRUE, data=NULL){
    if(is.null(data)){
        data <- getPlotDistributionData(gr, fds, type)
    }

    rac <- data$rocts + data$rcts
    rx <- x <- rac
    ry <- y <- data$rcts

    logpre <- ""
    logsuf <- ""
    if(plotLog){
        logpre <- "log10(1 + "
        logsuf <- ")"
        x <- log10(1 + x)
        y <- log10(1 + y)
    }

    heatscatter(x, y,
            main=getTitle("Heatscatter of raw counts", data$se, type),
            ylab=paste0(logpre, "raw counts of site of interest", logsuf),
            xlab=paste0(logpre, "raw all counts", logsuf)
    )

    # grid and diagonal
    abline(0,1,col="gray", lty="dotted")
    grid()


    a <- data$alpha
    b <- data$beta
    ab <- a + b
    fitx <- 0:as.integer(max(rac)*1.1)
    fity <- fitx * a/(a+b)
    fityvar <- (fitx * a * b * (fitx + ab))/(ab**2 * (ab + 1))
    fityvarbot <- fity - fityvar
    fityvartop <- fity + fityvar
    if(plotLog){
        fitx       <- log10(1+fitx)
        fity       <- log10(1+fity)
        fityvarbot <- log10(1+fityvarbot)
        fityvartop <- log10(1+fityvartop)
    }
    lines(fitx, fity, col="firebrick")
    lines(fitx, fityvartop, col="firebrick", lty="dotted")
    lines(fitx, fityvarbot, col="firebrick", lty="dotted")

    if(!is.null(sample)){
        names(x) <- samples(fds)
        names(y) <- samples(fds)
        sapply(sample, addSamplePoints, x=x, y=y)
    }

    if(plotLegend){
        legend("topleft", c("Model fit", "+/- Variance"),
                pch=20, lty=c(1,3), col="firebrick")
    }
}

#'
#' generic function to plot any value as a sample rank plot
#'
#' @noRd
plotSampleRank <- function(gr, fds, type, sample=NULL, delta=FALSE,
                    plotLog=FALSE, rmZeroCts=FALSE, data=NULL, ...){
    # get data
    if(is.null(data)){
        data <- getPlotDistributionData(gr, fds, type, rmZeroCts)
    }
    p <- data[[data$mapping[type]]]

    ylab=type
    if(delta){
        p <- p - median(p)
        ylab=paste0("delta_median( ", ylab, " )")
    }
    plot(sort(p),
        main=getTitle("Sample rank", data$se,
                paste0(ifelse(delta, "delta ", ""), type)),
        ylab=ylab, xlab="sample rank",
        pch=16, col="gray"
    )
    if(!is.null(sample)){
        sapply(sample, addSamplePoints, x=rank(p), y=p)
    }
}

#'
#' plot a given value against the counts
#'
plotValueVsCounts <- function(gr, fds, type, sample=NULL, delta=FALSE,
                    plotLog=FALSE, rmZeroCts=FALSE, data=NULL, ...){
    # get data
    if(is.null(data)){
        data    <- getPlotDistributionData(gr, fds, type, rmZeroCts)
    }
    p       <- data[[data$mapping[type]]]
    rcts    <- data$rcts
    rocts   <- data$rocts
    l10tcts <- log10(rcts + rocts)

    ylab=type
    if(delta){
        p <- data$deltaPsi
        ylab=paste0("delta_median( ", ylab, " )")
    }
    heatscatter(l10tcts, p,
            main=paste0("Heatscatter of ", type, " versus total raw counts"),
            xlab="log10(total raw counts)", ylab=type)

    if(!is.null(sample)){
        sapply(sample, addSamplePoints, x=l10tcts, y=p)
    }
}

#'
#' add the given sample as annotation to the plot
#'
#' @noRd
addSamplePoints <- function(x, y, sample, pch=20, col="red", ...){
    stopifnot(sample %in% names(x) && sample %in% names(y))

    points(x[sample], y[sample], pch=pch, col=col, ...)
    text(x[sample], y[sample], labels=sample, pos=3)
}

#'
#' create the title for all the plots
#'
#' @noRd
getTitle <- function(plotMainTxt, gr, psiType){
    paste0(plotMainTxt, " for type: ", psiType,
            "\nRange: ", seqnames(gr), ":", start(gr), "-", end(gr))
}

plotQQplot <- function(gr, fds, type, data=NULL){
    # get data
    if(is.null(data)){
        data    <- getPlotDistributionData(gr, fds, type)
    }
    o       <- sort(data$pvalues)
    e       <- ppoints(length(o))

    o[is.na(o)] <- 1
    plot(-log10(e), -log10(o), main="QQ-plot", pch=16)
    abline(0,1,col="firebrick")
    grid()
}

#' testing function
#' @noRd
testPlotting <- function(){
    testset <- FALSE
    if(!testset){
        fds  <- loadFraseRDataSet(
            "/s/project/fraser/analysis/datasets", "kremer-bader-et-al")
        grdt <- as.data.table(metadata(fds)[[3]])[p.adj < 1]
        gra   <- makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                grdt[order(p.adj)][
                        hgnc_symbol %in% c("MCOLN1", "TIMMDC1")][c(1,3,5)])
    } else {
        fds  <- annotateRanges(FraseR())
        grdt <- as.data.table(results(fds, fdrCut=1))
        gra   <- makeGRangesFromDataFrame(
            grdt[order(p.adj)][c(3,5,8)], keep.extra.columns = TRUE)
    }

    pdf("result-distribution-top50.pdf")
        idx2plot <- sample(1:min(3000, length(gra)), min(200, length(gra)))
        for(i in sort(idx2plot)){
            gr        <- gra[i]
            sampleIDs <- gra[i]$sampleID

            plotJunctionDistribution(fds=fds, gr=gr, sampleIDs=sampleIDs)
        }
    dev.off()
    browseURL("result-distribution-top50.pdf")

    debug(plotJunctionDistribution)
    debug(getPlotDistributionData)
    debug(plotSampleRank)
    debug(plotQQplot)
    plotJunctionDistribution(fds, curgr, valueVsCounts = TRUE, qqplot = TRUE)
}
