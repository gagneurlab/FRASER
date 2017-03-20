#'
#' Plot the results of the FraseR analysis pipeline
#' All three types are plotted into one HTML file
#' based on plotly.
#'
#' @export
#' @examples
#'      plotSampleResults(dataset, "sample1")
#'      plotSampleResults(dataset, "sample1", "result.html")
plotSampleResults <- function(dataset, sampleID=NULL,
                    file=NULL, browseIt=FALSE){
    require(FraseR)

    #
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")

    # if sample is empty create all plots for all samples
    if(is.null(sampleID)){
        samples2plot <- !is.na(sampleGroup(dataset@settings))
        sampleIDs2plot <- samples(dataset@settings[samples2plot])
        return(bplapply(sampleIDs2plot, plotSampleResults,
                dataset=dataset, file=file, BPPARAM=parallel(dataset@settings)
        ))
    }

    # check the rest
    stopifnot(sampleID %in% samples(dataset@settings))
    if(is.null(file)){
        outDir <- file.path(outputFolder(dataset@settings), "results")
        if(!is.null(outDir) && outDir != ""){
            file <- file.path(outDir, paste0(sampleID, "-FraseR-results.html"))
        }
    }

    # create folder if needed
    if(!is.null(file) && !dir.exists(dirname(file))){
        dir.create(dirname(file), recursive=TRUE)
    }

    # generate each sub plot
    psi3plot <- .plotVolcano(dataset, sampleID, "splitReads", "psi3", 1)
    psi5plot <- .plotVolcano(dataset, sampleID, "splitReads", "psi5", 2)
    psisplot <- .plotVolcano(dataset, sampleID, "nonSplicedReads", "sitePSI", 3)

    # combine plots
    mainplot <- plotly::subplot(psi3plot, psi5plot, psisplot,
            nrows = 3, shareX = FALSE, shareY = FALSE,
            titleX = TRUE, titleY = TRUE
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

    #
    if(!is.null(file)){
        saveWidget(mainplot, file=file)
        if(browseIt){
            browseURL(file)
        }
        return(file)
    }

    # return it
    return(mainplot)
}

#'
#' generate a volcano plot
#'
#' @noRd
.plotVolcano <- function(dataset, sampleID, readType, psiType, subID,
                    ylim=c(0,30), xlim=c(-5,5)){
    curSlot <- slot(dataset, readType)
    zscore  <- assays(curSlot)[[paste0("zscore_", psiType)]][,sampleID]
    pvalue  <- -log10(assays(curSlot)[[paste0("pvalue_", psiType)]][,sampleID])
    psival  <- assays(curSlot)[[psiType]][,sampleID]

    if("DelayedMatrix" %in% is(zscore)){
        zscore <- as.vector(zscore[,1])
    }
    if("DelayedMatrix" %in% is(pvalue)){
        pvalue <- as.vector(pvalue[,1])
    }
    if("DelayedMatrix" %in% is(psival)){
        psival <- as.vector(psival[,1])
    }

    # remove NAs from data
    toplot <- !is.na(zscore) & !is.na(pvalue) & !is.infinite(pvalue)

    # remove unsignificant data (to keep plotly responsive)
    # max 50k points
    xCutoff <- 1.5
    yCutoff <- 1.5
    unsigni <- abs(zscore[toplot]) < xCutoff & pvalue[toplot] < yCutoff
    if(sum(toplot[toplot] & !unsigni) < 50000){
        unsigni <- unsigni & rank(pvalue[toplot]) > 50000
    }
    toplot[toplot]  <- !unsigni


    # trim data to ylim and xlim
    pvalue2plot <- pvalue
    zscore2plot <- zscore
    pvalue2plot[toplot & pvalue2plot > max(ylim)] <- max(ylim)
    zscore2plot[toplot & zscore2plot > max(xlim)] <- max(xlim)
    zscore2plot[toplot & zscore2plot < min(xlim)] <- min(xlim)
    plotTraces <- list(
        "&#936; &#8804; 30%"=.na2false(psival[toplot]<=0.3),
        "30% < &#936; &#8804; 60%"=.na2false(psival[toplot]<=0.6 & psival[toplot]>0.3),
        "60% < &#936;"=.na2false(psival[toplot]>0.6),
        "&#936; &#8801; NA"=is.na(psival[toplot])
    )

    p <- plot_ly()
    for(i in seq_along(plotTraces)){
        t <- which(toplot)[plotTraces[[i]]]
        if(sum(t) == 0){
            next
        }
        p <- p %>% add_trace(type="scattergl",
            x=zscore2plot[t], y=pvalue2plot[t],
            mode = "markers",
            marker = list(
                color = pvalue2plot[t],
                cmin = 0,
                cmax = max(ylim),
                colorbar = list(
                    y = 0.8,
                    len = 0.4,
                    title = "-log<sub>10</sub> <i>P</i>-value"
                )
            ),
            #name = paste0(psiType, ": ", names(plotTraces)[i]),
            name = names(plotTraces)[i],
            legendgroup = names(plotTraces)[i],
            showlegend = psiType=="psi3",
            visible = ifelse(i<=2, TRUE, "legendonly"),
            text = paste0(
                "Symbol:     ", mcols(curSlot)$hgnc_symbol[t], "</br>",
                "Chromosome: ", seqnames(curSlot)[t],  "</br>",
                "Start:      ", start(curSlot)[t],     "</br>",
                "End:        ", end(curSlot)[t],       "</br>",
                "-log<sub>10</sub>(<i>P</i>-value): ", round(pvalue[t], 2), "</br>",
                "Z-score:    ", round(zscore[t], 2), "</br>",
                "PSI-value:  ", round(psival[t], 3)*100, "%</br>"
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

.na2false <- function(x){
    x[is.na(x)] <- FALSE
    return(x)
}

rerunPlot <- function(){
    source("./FraseR/R/plotResults.R")
    plotSampleResults(fds, "sample1")
}
