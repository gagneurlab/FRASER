#'
#' Plot the results of the FraseR analysis pipeline
#' All three types are plotted into one HTML file
#' based on plotly.
#'
#' @export
#' @examples 
#'      plotSampleResults(dataset, "sample1")
#'      plotSampleResults(dataset, "sample1", "result.html")
plotSampleResults <- function(dataset, sampleID, file=NULL){
    
    # generate each sub plot
    psi3plot    <- .plotVolcano(dataset, sampleID, "splitReads", "psi3", 1)
    psi5plot    <- .plotVolcano(dataset, sampleID, "splitReads", "psi5", 2)
    sitePSIplot <- .plotVolcano(dataset, sampleID, "nonSplicedReads", "sitePSI", 3)
    
    # combine plots
    mainplot <- plotly::subplot(psi3plot, psi5plot, sitePSIplot,
                  nrows = 3, shareX = T, shareY = F,
                  titleX = TRUE, titleY = TRUE
    ) %>% layout(showlegend = FALSE)
    
    #
    if(!is.null(file)){
        htmlwidgets::saveWidget(mainplot, file=file)
        return(NULL)
    } 
    
    # return it
    return(mainplot)
    
}

#'
#' generate a volcano plot
#'
#' @noRd
.plotVolcano <- function(dataset, sampleID, readType, psiType, subID, ylim=c(0,30), xlim=c(-5,5)){
    curSlot <- slot(dataset, readType)
    zscore  <- assays(curSlot)[[paste0("zscore_", psiType)]][,sampleID]
    pvalue  <- -log10(assays(curSlot)[[paste0("pvalue_", psiType)]][,sampleID])
    
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
    pvalue2plot[toplot & pvalue2plot > max(ylim)] <- max(ylim) - 1
    zscore2plot[toplot & zscore2plot > max(xlim)] <- max(xlim) - 1
    zscore2plot[toplot & zscore2plot < min(xlim)] <- min(xlim) + 1
    
    p <- plot_ly(x=zscore2plot[toplot], y=pvalue2plot[toplot], type="scatter",
                 mode = "markers",
                 marker = list(   
                     color = pvalue2plot[toplot],
                     cmin = 0,
                     cmax = max(ylim),
                     colorbar = list(yanchor = "middle")
                 ),
                 name = psiType,
                 text = paste0(
                     "Symbol:     ", mcols(curSlot)$hgnc_symbol[toplot], "</br>",
                     "Chromosome: ", seqnames(curSlot)[toplot],  "</br>",
                     "Start:      ", start(curSlot)[toplot],     "</br>",
                     "End:        ", end(curSlot)[toplot],       "</br>",
                     "-log<sub>10</sub>(<i>P</i>-value): ", round(pvalue[toplot], 2), "</br>",
                     "Z-score:    ", round(zscore[toplot], 2), "</br>"
                 )
    )
    p <- p %>% layout(
        xaxis=list(range=xlim, title="Z-score"),
        yaxis=list(range=ylim, title="-log<sub>10</sub><i>P</i>-value"),
        annotations=list(
            x=0, y=yCutoff/2, 
            text=paste("Removed points</br>",
                       round(sum(unsigni)/length(unsigni)*100, digits=2), 
                       "% of points"
            ),
            xref = paste0("x", subID),
            yref = paste0("y", subID),
            showarrow = FALSE
        )
    )
    
    # add removed area
    p <- p %>% layout(shapes = list(
            list(type = "rect",
                    fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
                    x0 = -xCutoff, x1 = xCutoff, xref = paste0("x", subID),
                    y0 = 0, y1 = yCutoff, yref = paste0("y", subID)
            )
    ))
    
    return(p)
}



