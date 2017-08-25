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
            rmZeroCts=FALSE, plotRank=paste0(c("", "delta_", "zscore_"), type),
            plotCounts=TRUE, plotValVsCounts=type, qqplot=TRUE,
            plotLegend=TRUE){
    stopifnot(is(fds, "FraseRDataSet"))
    if(is.data.table(gr)){
        gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE)
    }
    stopifnot(is(gr, "GRanges"))
    stopifnot(length(gr) == 1)

    if(is.null(sampleIDs) & "sampleID" %in% colnames(mcols(gr))){
        sampleIDs <-  gr$sampleID
    }

    # get number of plots
    numPlots <- sum(sapply(list(plotRank, plotValVsCounts),
                    function(x){ ifelse(is.logical(x), sum(x), length(x)) })) +
            sum(qqplot) + sum(plotCounts)

    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))

    data <- getPlotDistributionData(gr, fds, type, rmZeroCts)

    par(mfrow=c(ceiling(numPlots/3),3), cex=1)

    # plot sample rank if requested
    if(!(length(plotRank) == 0 | isFALSE(plotRank))){
        sapply(plotRank, plotSampleRank, gr=gr, fds=fds, sample=sampleIDs,
            rmZeroCts=rmZeroCts, data=data)
    }
    # plot counts
    if(isTRUE(plotCounts)){
        plotCountsAtSite(gr, fds, type=type, sample=sampleIDs, data=data,
                plotLegend=plotLegend)
    }
    # plot values versus counts
    if(!(length(valueVsCounts) == 0 | isFALSE(valueVsCounts))){
        sapply(plotValVsCounts, plotValueVsCounts, gr=gr, fds=fds, data=data,
                sampleIDs=sampleIDs, plotlog=TRUE, rmZeroCts=rmZeroCts)
    }
    # plot qq plot
    if(isTRUE(qqplot)){
        plotQQplot(gr, fds, type=type, data=data)
    }
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
        keepSample <- rcts + rocts > 0
        se <- se[,keepSample]
        rcts     = rcts[keepSample]
        rocts    = rocts[keepSample]
    }

    return(list(
        se       = seOfInterest,
        mapping  = mapping,
        pvalues  = as.matrix(assays(seOfInterest)[[pvaluename]])[1,],
        psi      = as.matrix(assays(seOfInterest)[[psiname]])[1,],
        deltaPsi = as.matrix(assays(seOfInterest)[[deltapsiname]])[1,],
        zscore   = as.matrix(assays(seOfInterest)[[zscorename]])[1,],
        rcts     = rcts,
        rocts    = rocts,
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
    # get data to plot
    if(is.null(data)){
        data <- getPlotDistributionData(gr, fds, type)
    }

    # raw counts
    rac <- data$rocts + data$rcts
    x <- rac
    y <- data$rcts
    transformFun <- function(x) x
    if(plotLog){
        transformFun <- function(x) suppressWarnings(log10(1+x))
    }

    # do we plot it in log?
    logpre <- logsuf <- ""
    if(plotLog){
        logpre <- "log10(1 + "
        logsuf <- ")"
    }
    xlab <- paste0(logpre, "raw all counts", logsuf)
    ylab <- paste0(logpre, "raw counts of site of interest", logsuf)

    if(type=="psiSite"){
        y <- data$rocts
        ylab <- paste0(logpre, "raw counts of spliced reads", logsuf)
    }

    # main heatscatter plot
    heatscatter(transformFun(x), transformFun(y), ylab=ylab, xlab=xlab,
            main=getTitle("Heatscatter of raw counts", data$se, type),
            xlim=c(0, max(transformFun(x))), ylim=c(0, max(transformFun(y))))

    # grid and diagonal
    abline(0,1,col="gray", lty="dotted")
    grid()

    # add prediction model
    a <- data$alpha
    b <- data$beta
    ab <- a + b
    fitx <- 0:as.integer(max(rac)*1.5)
    fity <- bbmean(fitx, a, b)
    lines(transformFun(fitx), transformFun(fity), col="firebrick")

    # add 50%, 25% and 10% lines
    curPlotPar <- data.table(fact=c(0.5, 0.25, 0.1), lty=c(2, 4, 3),
            name=c("PSI = 50%", "PSI = 25%", "PSI = 10%"))
    sapply(1:nrow(curPlotPar), function(idx){
        y2plot <- fitx * curPlotPar[idx, fact]
        lines(transformFun(fitx), transformFun(y2plot),
                lty=curPlotPar[idx, lty], col=adjustcolor("black", 0.7))
    })

    # add variance
    fityvar <- bbvariance(fitx, a, b)
    #scewness <- bbscewness(fitx, a, b)
    #scewness[1] <- 0
    #scewFactor <- ifelse(scewnewss < 0, 1+abs(scewnewss), 1/abs(scewnewss))
    sapply(c(-1, 1), function(varFactor) {
        y <- fity + varFactor * fityvar # * scewFactor
        lines(transformFun(fitx), transformFun(y),
                lty="dotted", col=adjustcolor("firebrick", 0.5))
    })

    # add sample annotation
    if(!is.null(sample)){
        names(x) <- samples(fds)
        names(y) <- samples(fds)
        sapply(sample, addSamplePoints, x=transformFun(x), y=transformFun(y))
    }

    # add legend if requested
    if(plotLegend){
        legend("topleft", c("Model fit", "+/- Variance", curPlotPar[,name]),
                pch=20, lty=c(1,3, curPlotPar[,lty]),
                col=c("firebrick",adjustcolor("firebrick", 0.5),
                        rep(adjustcolor("black", 0.7), 3)))
    }
}

#'
#' The implementation of the mean/variance/skewness of the
#' Beta-Binomial distribution
#'
#' https://en.wikipedia.org/wiki/Beta-binomial_distribution
#' @noRd
bbmean <- function(size, a, b){
    size * a / (a + b)
}

#'
#' The implementation of the mean/variance/skewness of the
#' Beta-Binomial distribution
#'
#' https://en.wikipedia.org/wiki/Beta-binomial_distribution
#' @noRd
bbvariance  <- function(size, a, b){
    numerator   <- size * a * b * (a + b + size)
    denuminator <- (a + b)^2 * (a + b + 1)
    return(numerator/denuminator)
}

#'
#' The implementation of the mean/variance/skewness of the
#' Beta-Binomial distribution
#'
#' https://en.wikipedia.org/wiki/Beta-binomial_distribution
#' @noRd
bbscewness <- function(size, a, b){
    presqrt  <- (a + b + 2 * size) * (b - a) / (a + b + 2)
    sqrtpart <- (1 + a + b) / (size * a * b * (size + a + b))
    return(presqrt * sqrt(sqrtpart))
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

plotQQplot <- function(gr, fds, type, data=NULL, maxOutlier=2){
    # get data
    if(is.null(data)){
        data <- getPlotDistributionData(gr, fds, type)
    }

    obs <- -log10(sort(data$pvalues))
    exp <- -log10(ppoints(length(obs)))
    len <- length(exp)

    obs[is.na(obs)] <- 0
    maxPoint <- max(c(exp, obs))
    ylim <- range(0, min(exp[1]*maxOutlier, maxPoint))

    # main plot area
    plot(NA, main="QQ-plot", xlim=range(exp), ylim=ylim,
            xlab=paste0(expression(log[10]), "(expected)"),
            ylab=paste0(expression(log[10]), "(observed)"))

    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    upper <- qbeta(0.025, 1:len, rev(1:len))
    lower <- qbeta(0.975, 1:len, rev(1:len))
    polygon(x=c(rev(exp), exp), y=-log10(c(rev(upper), lower)),
            col="gray", border="gray")

    # add the points
    points(exp, obs, pch=16)

    # diagonal and grid
    abline(0,1,col="firebrick")
    grid()

    # plot outliers
    outOfRange <- which(obs > max(ylim))
    if(length(outOfRange) > 0){
        points(exp[outOfRange], exp[1]*maxOutlier, pch=2, col='red')
    }

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
