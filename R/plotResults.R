#'
#'
#'@export
plotSampleResults <- function(dataset, sampleID){
    layout(matrix(1:4, nrow=2))
    .plotVolcano(dataset, sampleID, "splitReads", "psi3")
    .plotVolcano(dataset, sampleID, "splitReads", "psi5")
    .plotVolcano(dataset, sampleID, "nonSplicedReads", "sitePSI")
    plot(NA)
    
}

#'
#'
#' @noRd
.plotVolcano <- function(dataset, sampleID, readType, psiType){
    curSlot <- slot(dataset, readType)
    zscore  <- assays(curSlot)[[paste0("zscore_", psiType)]][,sampleID]
    pvalue  <- -log10(assays(curSlot)[[paste0("pvalue_", psiType)]][,sampleID])
    
    toplot <- !is.na(zscore) & !is.na(pvalue) & !is.infinite(pvalue) & pvalue > 2
    p <- plot_ly(x=zscore[toplot], y=pvalue[toplot], type="scattergl",
                 color = pvalue[toplot],
                 text = paste0(
                     "Symbol:     ", mcols(curSlot)$hgnc_symbol[toplot], "</br>",
                     "Chromosome: ", seqnames(curSlot)[toplot],  "</br>",
                     "Start:      ", start(curSlot)[toplot],     "</br>",
                     "End:        ", end(curSlot)[toplot],       "</br>"
                 )
    )
    
    htmlwidgets::saveWidget(as.widget(p), file="tmp.html")
    
    plot(zscore, -log10(pvalue), 
         ylab = "-log10(P-Value)", xlab = "Z-score",
         main = paste0("Volcano plot of ", toupper(psiType), "\n(Sample = ", sampleID, ")")
    )
}



