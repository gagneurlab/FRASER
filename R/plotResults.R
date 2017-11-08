#'
#' Plot the results of the FraseR analysis pipeline
#' All three types are plotted into one HTML file
#' based on plotly.
#'
#' @export
#' @examples
#'      fds <- FraseR()
#'      plotSampleResults(fds, "sample1")
#'      plotSampleResults(fds, "sample1", "result.html")
plotSampleResults <- function(fds, sampleID=NULL, file=NULL,
            dir=NULL, browseIt=FALSE){

    # if sample is empty create all plots for all samples
    if(is.null(sampleID)){
        if(!is.null(file))
            stop("Can't do multiple samples yet with file. File would be overriden!")

        samples2plot <- !is.na(condition(fds))
        sampleIDs2plot <- samples(fds)[samples2plot]
        return(bplapply(sampleIDs2plot, fds=fds, dir=dir, BPPARAM=parallel(fds),
            FUN=function(sampleID, fds, dir) {
                require(FraseR)
                plotSampleResults(fds, sampleID, dir=dir, browseIt=FALSE)
            }
        ))
    }

    # check the rest
    stopifnot(sampleID %in% samples(fds))
    if(is.null(file)){
        if(is.null(dir)){
            dir <- file.path(workingDir(fds), "results", nameNoSpace(fds))
        }
        file <- file.path(dir, paste0("FraseR-results-", sampleID, ".html"))
    }

    # create folder if needed
    if(!dir.exists(dirname(file))){
        dir.create(dirname(file), recursive=TRUE)
    }

    # create main plot
    mainplot <- createMainPlotFraseR(fds, sampleID)

    #
    saveWidget(mainplot[["plot"]], file=file)
    if(browseIt){
        browseURL(file)
    }
    return(file)
}

#'
#' generate the main FraseR plot
#' @noRd
createMainPlotFraseR <- function(fds, sampleID, source=NULL){
    ylim=c(0,30)
    xlim=c(-5,5)
    nxlim <- c(xlim[1]*1.05, xlim[2]*1.05)
    nylim <- c(ylim[1], ylim[2]*1.05)

    # generate each sub plot
    plotls <- lapply(c("psi3", "psi5", "psiSite"), function(x){
            plotVolcano(fds, sampleID, x, source=source, xlim=xlim, ylim=ylim)
    })

    # combine plots
    mainplot <- plotly::subplot(lapply(plotls, "[[", "plot"),
            nrows = 3, shareX = TRUE, shareY = TRUE
    )
    mainplot <- layout(mainplot,
        title = paste0("FraseR results for sample: <b>", sampleID, "</b>"),
        showlegend = TRUE,
        legend = list(
            x = 1,
            y = 0.2,
            title = "&#936; filter"
        ),
        yaxis =list(domain=c(0.70,1.00)),
        yaxis2=list(domain=c(0.35,0.65)),
        yaxis3=list(domain=c(0.00,0.30))
    )

    return(list(
        plot=mainplot,
        plotDF=rbindlist(lapply(plotls, "[[", "plotDF"))
    ))
}


#'
#' generate a volcano plot
#'
#' @noRd
plotVolcano <- function(fds, sampleID, psiType,
            ylim=c(0,30), xlim=c(-5,5), source=NULL){

    # extract values
    zscores  <- assays(fds)[[paste0("zscore_", psiType)]][,sampleID]
    pvalues  <- -log10(assays(fds)[[paste0("pvalue_", psiType)]][,sampleID])
    psivals  <- assays(fds)[[psiType]][,sampleID]
    counts   <- counts(fds, type=psiType, side="ofInterest")[,sampleID]
    ocounts  <- counts(fds, type=psiType, side="otherSide")[,sampleID]

    # convert from HDF5 to memory array
    zscores <- as.vector(zscores[,1])
    pvalues <- as.vector(pvalues[,1])
    psivals <- as.vector(psivals[,1])
    counts  <- as.vector(counts[,1])
    ocounts <- as.vector(ocounts[,1])

    # remove NAs from data
    toplot <- !is.na(zscores) & !is.na(pvalues) & !is.infinite(pvalues)

    # remove unsignificant data (to keep plotly responsive)
    # max 50k points
    xCutoff <- 1.5
    yCutoff <- 1.5
    unsigni <- abs(zscores[toplot]) < xCutoff & pvalues[toplot] < yCutoff
    toplot[toplot]  <- !unsigni

    # trim data to ylim and xlim
    pvalue2plot <- pvalues
    zscore2plot <- zscores
    pvalue2plot[toplot & pvalue2plot > max(ylim)] <- max(ylim)
    zscore2plot[toplot & zscore2plot > max(xlim)] <- max(xlim)
    zscore2plot[toplot & zscore2plot < min(xlim)] <- min(xlim)

    # filter bad betabinom results
    filteredRes <- na2false(
        (psivals <= 0.01 | psivals >= 0.99) &
        abs(zscores) <= xCutoff &
        counts + ocounts > 10000
    )[toplot]

    # traces to plot
    plotTraces <- list(
        "&#936; &#8804; 30%"       = na2false(psivals[toplot] <= 0.3),
        "30% < &#936; &#8804; 60%" = na2false(psivals[toplot] <= 0.6 & psivals[toplot] > 0.3),
        "60% < &#936;"             = na2false(psivals[toplot] > 0.6),
        "filtered"                 = filteredRes
    )

    plotDF <- data.table()
    p <- plot_ly(type="scatter", mode="markers", source=source)
    for(i in seq_along(plotTraces)){
        if(names(plotTraces)[i] != 'filtered'){
            t <- which(toplot)[plotTraces[[i]] & !filteredRes]
        } else {
            t <- which(toplot)[filteredRes]
        }
        if(length(t) == 0) {
            next
        }

        # generate data to plot
        tmpFds <- fds
        if(psiType == "psiSite") tmpFds <- nonSplicedReads(fds)
        if(is.null(mcols(fds, type=psiType)$hgnc_symbol)){
            mcols(fds, type=psiType)$hgnc_symbol <- NA
        }
        tmpData <- data.table(
            psiType   = psiType,
            traceNr   = i,
            traceName = names(plotTraces)[i],
            zscore2p  = zscore2plot,
            pvalue2p  = pvalue2plot,
            zscore    = zscores,
            pvalue    = pvalues,
            psivalue  = psivals,
            symbol    = mcols(fds, type=psiType)$hgnc_symbol,
            chr       = as.character(seqnames(tmpFds)),
            start     = start(tmpFds),
            end       = end(tmpFds),
            counts    = counts,
            ocounts   = ocounts
        )[t]
        tmpData[,pointNr:=1:.N]

        # show legend if we request shiny figure
        tmpShowlegend <- ifelse(is.null(source), psiType=="psiSite", TRUE)

        # save data in data.table object for later use
        plotDF <- rbind(plotDF, tmpData)

        # create trace
        p <- add_trace(p, data=tmpData, x=~zscore2p, y=~pvalue2p,
            marker = list(color = ~pvalue2p,
                cmin = 0, cmax = max(ylim),
                colorbar = list(y = 0.8, len = 0.4,
                        title = "-log<sub>10</sub>(<i>P</i>-value)"
                )
            ),
            name = names(plotTraces)[i],
            legendgroup = names(plotTraces)[i],
            showlegend = tmpShowlegend,
            visible = ifelse(i<=2, TRUE, "legendonly"),
            text = paste0(
                "Symbol:           ", tmpData$symbol,  "<br>",
                "Chromosome:       ", tmpData$chr,     "<br>",
                "Start:            ", tmpData$start,   "<br>",
                "End:              ", tmpData$end,     "<br>",
                "raw counts:       ", tmpData$counts,  "<br>",
                "raw other counts: ", tmpData$ocounts, "<br>",
                "-log<sub>10</sub>(<i>P</i>-value):   ", round(tmpData$pvalue, 2), "<br>",
                "Z-score:          ", round(tmpData$zscore, 2), "<br>",
                "PSI-value:        ", round(tmpData$psival, 3)*100, "%<br>"
            )
        )
    }

    nxlim <- c(xlim[1]*1.05, xlim[2]*1.05)
    nylim <- c(ylim[1], ylim[2]*1.05)

    # TODO: see: https://github.com/ropensci/plotly/issues/1019
    subID <- ""
    subY0axisAdjusted <- 0
    subY1axisAdjusted <- yCutoff
    if(is.null(source)){
        subID <- which(c("psi3", "psi5", "psiSite") == psiType)
        subY0axisAdjusted <- 6.0*(3-subID)
        subY1axisAdjusted <- yCutoff*3 + subY0axisAdjusted
    }

    p <- layout(p, showlegend = TRUE,
        xaxis=list(range=nxlim, title="Z-score"),
        # TODO P-value does not appear in italic
        yaxis=list(range=nylim, title=paste0("-log<sub>10</sub>(<i>P</i>-value)<br>", psiType)),
        shapes = list(list(
            type = "rect", fillcolor = "blue",
            line = list(color = "blue"), opacity = 0.3,
            x0 = -xCutoff, x1 = xCutoff,
            xref = paste0("x", subID), yref = paste0("y", subID),
            y0 = subY0axisAdjusted, y1 = subY1axisAdjusted
        ))
    )

    return(list(plot=p, plotDF=plotDF))
}

#'
#' plot count distribution
#'
plotCountsAtSite <- function(gr, fds, type, sample=NULL, plotLog=TRUE){
    se <- as(fds, "RangedSummarizedExperiment")
    if(type %in% "psiSite") {
        se <- nonSplicedReads(fds)
    }
    rcname <- paste0("rawCounts", toupper(checkReadType(fds, type)))
    rocname <- paste0("rawOtherCounts_",
        unlist(regmatches(type, gregexpr("psi(3|5|Site)", type, perl=TRUE)))
    )

    se2plot <- se[from(findOverlaps(granges(se), gr, type="equal"))]

    rcx  <- as.vector(assays(se2plot)[[rcname]])
    rocy <- as.vector(assays(se2plot)[[rocname]])
    logpre <- ""
    logsuf <- ""
    if(plotLog){
        logpre <- "log10(1 + "
        logsuf <- ")"
        rcx[rcx==0] <- 0.8
        rocy[rocy==0] <- 0.8
        rcx <- log10(rcx)
        rocy <- log10(rocy)
    }

    heatscatter(rcx, rocy,
            main=paste("Heatscatter of raw counts for type: ", type,
                    "\nRange: ", seqnames(se2plot), ":",
                    start(se2plot), "-", end(se2plot)
            ),
            xlab=paste0(logpre, "raw counts", logsuf),
            ylab=paste0(logpre, "raw other counts", logsuf)
    )
    if(!is.null(sample)){
        sapply(sample, addSamplePoints, fds=fds,
                x=rcx, y=rocy, pch=20, col="red"
        )
    }
}

addSamplePoints <- function(x,y,sample,fds, ...){
    idx <- which(samples(fds) == sample)
    points(x[idx], y[idx], ...)
    text(x[idx], y[idx], labels=sample, pos=3)
}
