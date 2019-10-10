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
            FUN=function(sID, fds, dir) {
                require(FraseR)
                plotSampleResults(fds, sID, dir=dir, browseIt=FALSE)
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
    plotls <- lapply(c("psi5", "psi3", "psiSite"), function(x){
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
#' Volcano plot
#'
#' Plots the p values over the delta psi values, known as volcano plot.
#' Visualizes per sample the outliers. By type and aggregate by
#' gene if requested.
#'
#' @export
plotVolcano <- function(fds, sampleID, type, minDeltaPsi=0.3, padjCutoff=0.05,
                    basePlot=TRUE, aggregate=FALSE,
                    main=paste0("Volcano plot: ", sampleID)){

    dt <- getPlottingDT(fds, axis="col", type=type, idx=sampleID,
            aggregate=aggregate)

    padj_line <- dt[padj < 0.05, min(5.5, -log10(pval))]
    if(length(padj_line) == 0){
        padj_line <- 5.5
    }

    dt[,aberrant:=FALSE]
    dt[abs(deltaPsi) >= minDeltaPsi & padj <= padjCutoff, aberrant:=TRUE]

    g <- ggplot(dt, aes(x=deltaPsi, y=-log10(pval), color=aberrant, text=paste0(
                "SampleID: ", sampleID, "<br>",
                "featureID: ", featureID, "<br>",
                "Raw count (K): ", k, "<br>",
                "Raw total count (N): ", n, "<br>",
                "p value: ", signif(pval, 5), "<br>",
                "delta Psi: ", round(deltaPsi, 2), "<br>",
                "Type: ", type))) +
        geom_point(aes(alpha=ifelse(isTRUE(aberrant), 1, 0.5))) +
        xlab(expression(paste(Delta, Psi))) +
        ylab(expression(paste(-log[10], "(p value)"))) +
        geom_vline(xintercept=c(-minDeltaPsi, minDeltaPsi),
                color="firebrick", linetype=2) +
        geom_hline(yintercept=padj_line, color="firebrick", linetype=4) +
        ggtitle(main) +
        theme_bw() +
        theme(legend.position="none") +
        scale_color_manual(values=c("gray70", "firebrick"))

    if(isFALSE(basePlot)){
        g <- g + xlab("delta Psi") +
            ylab("-log[10](p value)")
        return(ggplotly(g, tooltip="text"))
    }
    g
}


#'
#' Number of aberrant events per sample
#'
#' Plot the number of aberrant events per samples
#'
#' @export
plotAberrantPerSample <- function(fds, main, padjCutoff=0.05, zScoreCutoff=NA,
                    deltaPsiCutoff=0.1, col="#1B9E77", type=c("psi3", "psi5", "psiSite"),
                    yadjust=c(1.2, 1.2), labLine=c(3.5, 3), ymax=NULL,
                    ylab="#Aberrantly spliced events", labCex=par()$cex, ...){

    type <- match.arg(type)

    if(missing(main)){
        main <- paste('Aberrant events per sample (', type, ')')
    }

    count_vector <- sort(aberrant(fds, by="sample", type=type,
            padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
            deltaPsiCutoff=deltaPsiCutoff, ...))

    ylim <- c(0.4, max(1, count_vector)*1.1)
    if(!is.null(ymax)){
        ylim[2] <- ymax
    }
    replace_zero_unknown <- 0.5
    ticks <- c(replace_zero_unknown, signif(10^seq(
        from=0, to=round(log10(max(1, count_vector))), by=1/3), 1))

    labels_for_ticks <- sub(replace_zero_unknown, '0', as.character(ticks))

    bp <- barplot2(
        replace(count_vector, count_vector==0, replace_zero_unknown),
        log='y', ylim=ylim, names.arg='', xlab='', plot.grid=TRUE,
        grid.col='lightgray', ylab='', yaxt='n', border=NA, xpd=TRUE,
        col=col, main=main)

    n_names <- floor(length(count_vector)/20)
    xnames= seq_len(n_names*20)
    axis(side=1, at= c(0,bp[xnames,]), labels= c(0,xnames))
    axis(side=2, at=ticks, labels= labels_for_ticks, ylog=TRUE, las=2)

    # labels
    mtext('Sample rank', side=1, line=labLine[1], cex=labCex)
    mtext(ylab, side=2, line=labLine[2], cex=labCex)

    # legend and lines
    hlines <- c(Median=max(replace_zero_unknown, median(count_vector)),
            Quantile90=quantile(count_vector, 0.9, names=FALSE))
    color_hline <- c('black','black')
    abline(h=hlines, col=color_hline)
    text(x=c(1,1), y=hlines*yadjust, col=color_hline, adj=0,
         labels=c('Median', expression(90^th ~ 'percentile')))

    box()
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
