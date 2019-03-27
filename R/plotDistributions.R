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
#' @examples
#'
#' fds <- createTestFraseRDataSet()
#'
#' plotJunctionDistribution(fds, results(fds)[1])
#' plotJunctionDistribution(fds, results(fds)[3])
#'
#' @export
plotJunctionDistribution <- function(fds, gr, type=gr$type, sampleIDs=NULL,
            rmZeroCts=FALSE, plotRank=paste0(c("", "delta_", "zscore_"), type),
            plotCounts=TRUE, plotValVsCounts=type, qqplot=TRUE,
            plotLegend=TRUE, cex=1, mfrow=3, ...){
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

    data <- getPlotDistributionData(gr=gr, fds=fds, type=type, rmZeroCts)

    par(mfrow=c(ceiling(numPlots/mfrow),ifelse(numPlots < mfrow, numPlots, mfrow)), cex=cex)

    # plot sample rank if requested
    if(!(length(plotRank) == 0 | isFALSE(plotRank))){
        for(i in seq_along(plotRank)){
            plotSampleRank(gr=gr, fds=fds, type=plotRank[i],
                    sampleIDs=sampleIDs, rmZeroCts=rmZeroCts, data=data)
        }
    }
    # plot counts
    if(isTRUE(plotCounts)){
        plotCountsAtSite(gr, fds, type=type, sampleIDs=sampleIDs, data=data,
                plotLegend=plotLegend)
    }
    # plot values versus counts
    if(!(length(plotValVsCounts) == 0 | isFALSE(plotValVsCounts))){
        sapply(plotValVsCounts, plotValueVsCounts, gr=gr, fds=fds, data=NULL,
               sampleIDs=sampleIDs, plotlog=TRUE, rmZeroCts=FALSE)
    }
    # plot qq plot
    if(isTRUE(qqplot)){
        plotQQplot(gr, fds, type=type, data=data, ...)
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

    alpha <- beta <- NA
    if(paste0(type, "_alpha") %in% colnames(mcols(seOfInterest))){
        alpha <- mcols(seOfInterest)[,paste0(type, "_alpha")]
        beta  <- mcols(seOfInterest)[,paste0(type, "_beta")]
    }

    return(list(
        se       = seOfInterest,
        mapping  = mapping,
        rcts     = rcts,
        rocts    = rocts,
        pvalues  = as.matrix(null2na(assays(seOfInterest)[[pvaluename]]))[1,],
        psi      = as.matrix(null2na(assays(seOfInterest)[[psiname]]))[1,],
        deltaPsi = as.matrix(null2na(assays(seOfInterest)[[deltapsiname]]))[1,],
        zscore   = as.matrix(null2na(assays(seOfInterest)[[zscorename]]))[1,],
        alpha    = alpha,
        beta     = beta
    ))
}

#'
#' This returns the zero value replacement in log plots for counts
#'
#' @noRd
getZeroReplacement <- function(..., currZeroVal){
    if(!is.na(currZeroVal) && numeric(currZeroVal) &&
                currZeroVal < 1 && currZeroVal > 0.2){
        return(currZeroVal)
    }
    maxVal <- max(unlist(...))
    if(maxVal < 1000){
        return(0.64)
    }
    if(maxVal < 5000){
        return(0.68)
    }
    return(0.7)
}

#'
#' plot count distribution
#'
#' @noRd
plotCountsAtSite <- function(gr, fds, type, sampleIDs=NULL, plotLegend=TRUE,
                    data=NULL, zeroVal=NA){
    # get data to plot
    if(is.null(data)){
        data <- getPlotDistributionData(gr, fds, type)
    }
    if(isTRUE(plotLegend)){
        plotLegend <- "topleft"
    }

    # raw counts
    x <- data$rocts + data$rcts
    y <- data$rcts
    zeroVal <- getZeroReplacement(x, y, currZeroVal=zeroVal)
    x[x == 0] <- zeroVal
    y[y == 0] <- zeroVal

    # x and y lables
    xlab <- "All read counts"
    ylab <- "Read counts for side of interest"
    main <- getTitle("Heatscatter of raw counts", data$se, type)

    if(type=="psiSite"){
        y <- x - y
        ylab <- "Non spliced read counts"
    }

    # main heatscatter plot
    heatscatter(x, y, ylab=ylab, xlab=xlab, main=main, log="xy")

    # grid and diagonal
    abline(0,1,col="gray", lty="dotted")
    grid()

    if(any(x < 1)){
        abline(v=zeroVal, col="gray", lty="dotted")
        axis(side=1, at=zeroVal, labels = "0", tick = TRUE)
    }
    if(any(y < 1)){
        abline(h=zeroVal, col="gray", lty="dotted")
        axis(side=2, at=zeroVal, labels = "0", tick = TRUE)
    }

    # add prediction model
    a <- data$alpha
    b <- data$beta
    ab <- a + b
    fitx <- 1:as.integer(max(x)*1.5)
    fity <- bbmean(fitx, a, b)
    lines(c(zeroVal, fitx), c(zeroVal, fity), col="firebrick")

    # add 50%, 25% and 10% lines
    curPlotPar <- data.table(fact=c(0.5, 0.25, 0.1), lty=c(2, 4, 3),
            name=c("PSI = 50%", "PSI = 25%", "PSI = 10%"))
    sapply(1:nrow(curPlotPar), function(idx){
        y2plot <- fitx * curPlotPar[idx, fact]
        lines(fitx, y2plot, lty=curPlotPar[idx, lty],
                col=adjustcolor("black", 0.7))
    })

    # add variance
    fityvar <- bbvariance(fitx, a, b)
    sapply(c(-1, 1), function(varFactor) {
        y <- fity + varFactor * fityvar # * scewFactor
        lines(fitx, y, lty="dotted", col=adjustcolor("firebrick", 0.5))
    })

    # add sample annotation
    if(!is.null(sampleIDs)){
        names(x) <- samples(fds)
        names(y) <- samples(fds)
        sapply(sampleIDs, addSamplePoints, x=x, y=y)
    }

    # add legend if requested
    if(isScalarCharacter(plotLegend)){
        legend(plotLegend, c("Model fit", "+/- Variance", curPlotPar[,name]),
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
plotSampleRank <- function(gr, fds, type, sampleIDs=NULL, delta=FALSE,
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
    if(!is.null(sampleIDs)){
        sapply(sampleIDs, addSamplePoints, x=rank(p), y=p)
    }
}

#'
#' plot a given value against the counts
#' @noRd
plotValueVsCounts <- function(gr, fds, type, sampleIDs=NULL, delta=FALSE,
                    main=paste0(type, " versus total raw counts"),
                    xlab="Number of all observed split reads", ylab=NULL,
                    zeroVal=NA, rmZeroCts=FALSE, data=NULL, ...){
    # get data
    if(is.null(data)){
        data <- getPlotDistributionData(gr, fds, type, rmZeroCts)
    }
    val <- data[[data$mapping[type]]]
    tcts <- data$rcts + data$rocts
    zeroVal <- getZeroReplacement(tcts, currZeroVal=zeroVal)
    tcts[tcts == 0] <- zeroVal

    if(is.null(ylab)){
        ylab=type
        if(delta){
            val <- data$deltaPsi
            ylab=paste0("delta_median( ", ylab, " )")
        }
    }

    heatscatter(tcts, val, main="", xlab=xlab, log="x", ylab=ylab)
    title(main=main)
    grid(equilogs = FALSE)

    if(any(tcts < 1)){
        abline(v=zeroVal, col="gray", lty="dotted")
        axis(side=1, at=zeroVal, labels = "0", tick = TRUE)
    }
    if(!is.null(sampleIDs)){
        sapply(sampleIDs, addSamplePoints, x=tcts, y=val)
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

#'
#' qqplot
#'
#' @noRd
plotQQplot <- function(gr=NULL, fds=NULL, type=NULL, data=NULL, maxOutlier=2,
                    conf.alpha=0.05, sample=FALSE, breakTies=TRUE,
                    legendPos="bottomright"){
    if(isScalarLogical(conf.alpha)){
        conf.alpha <- ifelse(isTRUE(conf.alpha), 0.05, NA)
    }

    # get data
    if(is.null(data)){
        if(is.null(gr) | is.null(fds) | is.null(type)){
            stop("If data is not provided gr, fds and type needs to be passed on.")
        }
        data <- getPlotDistributionData(gr, fds, type)
    }

    # points
    obs <- -log10(sort(data$pvalues))
    obs[is.na(obs)] <- 0
    if(length(unique(obs)) < 2 | length(obs) < 2){
        warning("There are no pvalues or all are NA!")
        return(FALSE)
    }
    if(breakTies){
        obs <- breakTies(obs, logBase=10, decreasing=TRUE)
    }
    exp <- -log10(ppoints(length(obs)))
    len <- length(exp)

    # limits for nice plotting
    maxPoint <- max(c(exp, obs))
    ylim <- range(0, min(exp[1]*maxOutlier, maxPoint))

    # main plot area
    plot(NA, main=mainName, xlim=range(exp), ylim=ylim,
         xlab=expression(-log[10] ~  "(expected P-value)"),
         ylab=expression(-log[10] ~ "(observed P-value)"))


    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(is.numeric(conf.alpha)){
        getY <- function(x, exp){
            x1 <- exp[2]
            x2 <- exp[1]
            y1 <- -log10(x[2])
            y2 <- -log10(x[1])
            m <- (y2-y1)/(x2-x1)
            return(10^-(y1 + m*((x2+1)-x1)))
        }
        upper <- qbeta(conf.alpha/2,   1:len, rev(1:len))
        lower <- qbeta(1-conf.alpha/2, 1:len, rev(1:len))
        polygon(col="gray", border="gray", x=c(rev(exp), max(exp)+c(1,1), exp),
                y=-log10(c(
                        rev(upper), getY(upper, exp), getY(lower, exp), lower)))
    }

    # grid
    grid()

    plotPoint <- TRUE
    if(isTRUE(sample)){
        lo <- length(obs)
        plotPoint <- 1:lo %in% unique(c(1:min(lo, 100), sort(sample(1:lo,
                size = min(lo, 30000), prob=log10(1+lo:1)/sum(log10(1+lo:1))))))
    }

    # add the points
    outOfRange <- obs > max(ylim)
    points(exp[plotPoint & !outOfRange], obs[plotPoint & !outOfRange], pch=16)

    # diagonal and grid
    abline(0,1,col="firebrick")

    # plot outliers
    if(sum(outOfRange) > 0){
        points(exp[plotPoint & outOfRange], rep(max(ylim), sum(outOfRange)),
                pch=2, col='red')
    }

    if(is.numeric(conf.alpha)){
        legend(legendPos, paste0("CI (alpha = ",
                signif(conf.alpha, 2), ")"), lty=1, lwd=7, col="gray")
    }
}

#' sample qq
#'
#' @noRd
plotSampleQQ <- function(fds, type=c("psi5", "psi3", "psiSite"), sample=TRUE,
                    ...){
    pvals <- sapply(type, function(x){
        readType <- whichReadType(fds, x)
        tested <- na2false(mcols(fds, type=x)[,paste0(x, "_tested")])
        as(assays(fds[tested,by=readType])[[paste0('pvalue_', x)]], "matrix")
    })
    plotQQplot(data=list(pvalues=as.vector(unlist(pvals))), sample=sample, ...)
}

#'
#' breaks ties in a qq plot to get a better distributed p-value plot
#' @noRd
breakTies <- function(x, logBase=10, decreasing=TRUE){
    intervals <- sort(unique(c(0, x)))
    idxintervals <- findInterval(x, intervals)
    for(idx in as.integer(names(which(table(idxintervals) > 1)))){
        if(is.numeric(logBase)){
            minval <- logBase^-intervals[idx+1]
            maxval <- logBase^-intervals[idx]
            rand   <- runif(sum(idxintervals==idx), minval, maxval)
            rand   <- -log(rand, logBase)
        } else {
            minval <- intervals[idx]
            maxval <- intervals[idx+1]
            rand   <- runif(sum(idxintervals==idx), minval, maxval)
        }
        x[idxintervals==idx] <- rand
    }
    if(!is.na(decreasing)){
        x <- sort(x, decreasing=TRUE)
    }
}

#' testing function
#' @noRd
testPlotting <- function(){
    testset <- FALSE
    if(!testset){
        fds  <- loadFraseRDataSet(
            "/s/project/fraser/analysis/datasets", "kremer-bader-et-al")
        grdt <- as.data.table(metadata(fds)[[6]])[p.adj < 1 & abs(deltaPsi) > 0.1]
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
    debug(plotCountsAtSite)
    debug(plotValueVsCounts)
    debug(plotSampleRank)
    debug(plotQQplot)
    plotJunctionDistribution(fds, curgr, valueVsCounts = TRUE, qqplot = TRUE)
    plot(1,1)
    plotJunctionDistribution(fds, gra[1], rmZeroCts = FALSE)
    plotCountsAtSite(gra[1], fds, gra[1]$type)
    plotValueVsCounts(gra[1], fds, gra[1]$type, plotLog=TRUE, rmZeroCts=FALSE)
}

#'
#' Plot count correlation
#'
#' @export
plotCountCorHeatmap <- function(fds, type=c("psi5", "psi3", "psiSite"),
            logit=FALSE, topN=50000, minMedian=1, main=NULL,
            show_rownames=FALSE, show_colnames=FALSE,
            annotation_col=NA, annotation_row=NA, ...){

    type <- match.arg(type)

    kmat <- as.matrix(counts(fds, type=type, side="ofIn"))
    nmat <- kmat + as.matrix(counts(fds, type=type, side="other"))

    expRowsMedian <- rowMedians(kmat) >= minMedian
    expRowsMax    <- rowMax(kmat)     >= 10
    table(expRowsMax & expRowsMedian)

    skmat <- kmat[expRowsMax & expRowsMedian,]
    snmat <- nmat[expRowsMax & expRowsMedian,]

    xmat <- (skmat + 1)/(snmat + 2)
    if(isTRUE(logit)){
        xmat <- qlogis(xmat)
    }
    xmat_rc    <- xmat - rowMeans(xmat)
    xmat_rc_sd <- rowSds(xmat_rc)
    cormat <- cor(xmat_rc[rank(xmat_rc_sd) >= length(xmat_rc_sd) - topN,])

    if(is.character(annotation_col)){
        annotation_col <- getColDataAsDFFactors(fds, annotation_col)
    }
    if(is.character(annotation_row)){
        annotation_row <- getColDataAsDFFactors(fds, annotation_row)
    }

    if(is.null(main)){
        if(isTRUE(logit)){
            main <- paste0("Row-centered Logit(PSI) correlation (", type, ")")
        } else {
            main <- paste0("Row-centered PSI correlation (", type, ")")
        }
    }

    pheatmap(cormat, show_rownames=show_rownames, show_colnames=show_colnames,
             main=main, annotation_col=annotation_col,
             annotation_row=annotation_row, ...,
             breaks=seq(-1, 1, length.out=50),
             color=colorRampPalette(colors=rev(brewer.pal(11, "RdBu")))(50)
    )
}


getColDataAsDFFactors <- function(fds, names){
    tmpDF <- data.frame(colData(fds)[,names])
    colnames(tmpDF) <- names
    for(i in names){
        if(any(is.na(tmpDF[, i]))){
            tmpDF[,i] <- as.factor(paste0("", tmpDF[,i]))
        }
    }
    rownames(tmpDF) <- rownames(colData(fds))
    return(tmpDF)
}
