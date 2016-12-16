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
    zscore <- assays(slot(dataset, readType))[[paste0("zscore_", psiType)]][,sampleID]
    pvalue <- assays(slot(dataset, readType))[[paste0("pvalue_", psiType)]][,sampleID]
    plot(zscore, -log10(pvalue), 
         ylab = "-log10(P-Value)", xlab = "Z-score",
         main = paste0("Volcano plot of ", toupper(psiType), "\n(Sample = ", sampleID, ")")
    )
}



