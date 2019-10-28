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
createMainPlotFraseR <- function(fds, sampleID, ...){
    ylim=c(0,30)
    xlim=c(-5,5)
    nxlim <- c(xlim[1]*1.05, xlim[2]*1.05)
    nylim <- c(ylim[1], ylim[2]*1.05)

    # generate each sub plot
    plotls <- lapply(c("psi5", "psi3", "psiSite"), function(x){
            plotVolcano(fds, sampleID, type=x, ...)
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
