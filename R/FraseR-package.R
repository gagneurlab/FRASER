#'
#' FraseR: A package providing a workflow to detect aberrant splicing events
#'     in RNA sequencing data in a rare disease cohort
#'
### PACKAGE NAME
#'
#' @name FraseR
#'
#'
### Parallel computing
#'
#' @importFrom parallel mclapply
#' @import BiocParallel
#' @importFrom pcaMethods pca loadings
#'
### GRange/Experiment/bamFile packages
#'
#' @importFrom BiocGenerics updateObject counts
#' @import GenomicAlignments
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps from to IRanges
#' @importFrom Rsubread featureCounts
#' @importFrom Rsamtools ScanBamParam scanBamHeader bamMapqFilter
#'          bamWhich bamWhich<-
#'
#'
### Annotation
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#'
### Plotting
#'
#' @importFrom graphics plot par
#' @importFrom plotly plot_ly subplot layout add_trace
#' @importFrom gplots col2hex
#' @importFrom htmlwidgets saveWidget
#' @importFrom LSD heatscatter
#' @importFrom pheatmap pheatmap
#'
### Shiny App
#'
## @import shiny
#' @importFrom shiny column fluidPage fluidRow h3 h1 htmlOutput mainPanel
#'          navbarPage numericInput p plotOutput renderText selectInput
#'          sidebarLayout sidebarPanel shinyApp tabPanel textInput titlePanel
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly plotlyOutput renderPlotly
#'
### Data handling
#'
#' @importFrom data.table data.table as.data.table is.data.table := fread
#'          setnames
#' @import tidyr
#' @importFrom HDF5Array writeHDF5Array path
#' @importFrom rhdf5 h5ls
#'
### P-Value calculation
#'
#' @importFrom stats sd rbinom fisher.test na.omit p.adjust ppoints qbeta rnorm
#' @importFrom VGAM rbetabinom vglm Coef pbetabinom pbetabinom.ab betabinomial
#'          dbetabinom.ab dbetabinom
#'
### Miscelenious functions
#'
#' @importFrom BBmisc isScalarCharacter isScalarLogical chunk %nin%
#'          isScalarInteger isFALSE is.error
#' @importFrom R.utils renameFile
#' @importFrom tools file_path_as_absolute
#' @importFrom methods as callNextMethod is new slot slot<- validObject
#' @importFrom utils browseURL capture.output sessionInfo
#'
#'
#
### To be added into the functions above
#
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom grDevices dev.off adjustcolor pdf
#' @importFrom graphics abline axis grid legend lines title text points polygon
#' @importFrom plotly event_data
#' @importFrom GenomeInfoDb seqlevels<- seqlevels seqlengths
#'          keepStandardChromosomes
#' @importFrom shiny renderUI reactive renderPrint renderPlot
#' @importFrom DelayedArray rowMaxs rowMeans path<-
#' @importFrom data.table rbindlist
#' @importFrom DelayedMatrixStats rowMedians rowSds
#' @importFrom stats runif median quantile
#' @importFrom extraDistr rdirmnom
#' @importFrom PRROC pr.curve
#'
#' @useDynLib FraseR
#'



NULL

#' TODO: Should be removed in the end!
#' @noRd
globalVariables(c(".N", ".asDataFrame", "End", "FN", "HTML", "Start", "TP",
        "deltaPsi", "curgr", "gene", "lty", "hgnc_symbol", "id",
        "ldat", "p.adj", "pval", "pvalue", "shinyFds", "shinyFdsRes",
        "sampleID", "sampleGroup", "chr", "symbol", "type", "pseudocount"),
        package="FraseR")


options("FraseR.pseudoCount"=1)


