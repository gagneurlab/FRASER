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
plotSampleResults <- function(fds, sampleID=NULL, file=NULL, browseIt=FALSE){

    # check input
    stopifnot(class(fds) == "FraseRDataSet")

    # if sample is empty create all plots for all samples
    if(is.null(sampleID)){
        if(!is.null(file))
            stop("Can't do multiple samples yet with file. File would be overriden!")
        samples2plot <- !is.na(condition(fds))
        sampleIDs2plot <- samples(fds)[samples2plot]
        return(bplapply(sampleIDs2plot, plotSampleResults,
                        fds=fds, BPPARAM=parallel(fds)
        ))
    }

    # check the rest
    stopifnot(sampleID %in% samples(fds))
    if(is.null(file)){
        outDir <- file.path(workingDir(fds), "results", nameNoSpace(fds))
        file <- file.path(outDir, paste0(sampleID, "-FraseR-results.html"))
    }

    # create folder if needed
    if(!is.null(file) && !dir.exists(dirname(file))){
        dir.create(dirname(file), recursive=TRUE)
    }

    # create main plot
    mainplot <- createMainPlotFraseR(fds, sampleID)

    #
    saveWidget(mainplot, file=file)
    if(browseIt){
        browseURL(file)
    }
    return(file)
}

#'
#' generate the main FraseR plot
#' @noRd
createMainPlotFraseR <- function(fds, sampleID){

    # generate each sub plot
    psi3plot <- plotVolcano(fds, sampleID, "psi3")
    psi5plot <- plotVolcano(fds, sampleID, "psi5")
    psisplot <- plotVolcano(fds, sampleID, "psiSite")

    # combine plots
    mainplot <- subplot(psi3plot, psi5plot, psisplot, nrows = 3,
            shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE
    ) %>% layout(showlegend = TRUE,
            #xaxis=list(title="Z-score"),
            # TODO P-value does not appear in italic
            #yaxis=list(title="-log10 P-value"),

            legend = list(
                    x = 1,
                    y = 0.1,
                    title = "&#936; filter"
            ),
            title = paste0("FraseR results for sample: <b>", sampleID, "</b>")
    )

    return(mainplot)
}


#'
#' generate a volcano plot
#'
#' @noRd
plotVolcano <- function(fds, sampleID, psiType, ylim=c(0,30), xlim=c(-5,5)){
    # extract values
    zscores  <- assays(fds)[[paste0("zscore_", psiType)]][,sampleID]
    pvalues  <- -log10(assays(fds)[[paste0("pvalue_", psiType)]][,sampleID])
    psivals  <- assays(fds)[[psiType]][,sampleID]

    # convert from HDF5 to memory array
    zscores <- as.vector(zscores[,1])
    pvalues <- as.vector(pvalues[,1])
    psivals <- as.vector(psivals[,1])

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

    # traces to plot
    plotTraces <- list(
        "&#936; &#8804; 30%" =       na2false(psivals[toplot]<=0.3),
        "30% < &#936; &#8804; 60%" = na2false(psivals[toplot]<=0.6 & psivals[toplot]>0.3),
        "60% < &#936;" =             na2false(psivals[toplot]>0.6),
        "&#936; &#8801; NA" =        is.na(psivals[toplot]) # currently emtpy
    )

    p <- plot_ly()
    for(i in seq_along(plotTraces)){
        t <- which(toplot)[plotTraces[[i]]]
        if(length(t) == 0){
            next
        }

        # generate data to plot
        tmpFds <- fds
        if(psiType == "psiSite") tmpFds <- nonSplicedReads(fds)
        if(is.null(mcols(fds, type=psiType)$hgnc_symbol)){
            mcols(fds, type=psiType)$hgnc_symbol <- NA
        }
        tmpData <- data.table(
            zscore2p = zscore2plot,
            pvalue2p = pvalue2plot,
            zscore   = zscores,
            pvalue   = pvalues,
            psivalue = psivals,
            symbol   = mcols(fds, type=psiType)$hgnc_symbol,
            chr      = as.character(seqnames(tmpFds)),
            start    = start(tmpFds),
            end      = end(tmpFds),
            counts   = as.vector(counts(
                    fds, type=psiType, side="ofInterest")[,sampleID]),
            ocounts  = as.vector(counts(
                    fds, type=psiType, side="otherSide")[,sampleID])
        )[t]

        p <- p %>% add_trace(data=tmpData, type="scattergl", mode = "markers",
            x=~zscore2p, y=~pvalue2p,
            marker = list(color = ~pvalue2p,
                    cmin = 0, cmax = max(ylim),
                    colorbar = list(y = 0.8, len = 0.4,
                            title = "-log<sub>10</sub> <i>P</i>-value"
                    )
            ),
            name = names(plotTraces)[i],
            legendgroup = names(plotTraces)[i],
            showlegend = psiType=="psi3",
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

    p <- p %>% layout(
        xaxis=list(range=nxlim, title=ifelse(psiType=="sitePSI", "Z-score", "")),
        # TODO P-value does not appear in italic
        yaxis=list(range=nylim, title=paste0("-log10 P-value for ", psiType)),
        showlegend = FALSE
        #annotations=list(
        #    x=0,
        #    y=yCutoff/2,
        #    text="",
        #text=paste("Removed points</br>",
        #    round(sum(unsigni)/length(unsigni)*100, digits=2),
        #    "% of points"
        #),
        #    xref = paste0("x", subID),
        #    yref = paste0("y", subID),
        #    showarrow = FALSE
        #)
    )

    # add removed area
    #p <- p %>% layout(shapes = list(
    #        list(type = "rect",
    #                fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
    #                x0 = -xCutoff, x1 = xCutoff, xref = paste0("x", subID),
    #                y0 = 0, y1 = yCutoff, yref = paste0("y", subID)
    #        )
    #))

    return(p)
}

rerunPlot <- function(){
    library(htmlwidgets)
    library(plotly)
    library(FraseR)

    fds
    sampleID <- "SRR1087369"
    psiType <- "psi3"
    ylim=c(0,30)
    xlim=c(-5,5)

    na2false <- FraseR:::na2false
    nameNoSpace <- FraseR:::nameNoSpace

    plotSampleResults(fds, sampleID)

    source("./FraseR/R/plotResults.R")
    createMainPlotFraseR(fds, "sample1")
}
