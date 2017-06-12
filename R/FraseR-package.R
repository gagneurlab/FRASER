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
#' @importFrom BiocParallel SerialParam MulticoreParam
#'
#'
### GRange/Experiment/bamFile packages
#'
#' @import GenomicAlignments
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps from to IRanges
#' @importFrom Rsubread featureCounts
#' @importFrom Rsamtools ScanBamParam scanBamHeader
#' @importFrom Rsamtools bamMapqFilter bamWhich bamWhich<-
#'
#'
### Annotation
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#'
### Plotting
#'
#' @importFrom graphics plot
#' @importFrom plotly plot_ly subplot layout add_trace
#' @importFrom htmlwidgets saveWidget
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom DelayedArray rowMeans
#' @importFrom rhdf5 h5ls
#'
### Shiny App
#'
#' @importFrom shiny fluidPage titlePanel h1 fluidRow column numericInput h3
#' @importFrom shiny htmlOutput textInput sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny tabPanel navbarPage shinyApp selectInput p
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly plotlyOutput renderPlotly
#'
### Data handling
#'
#' @importFrom data.table data.table as.data.table is.data.table
#' @importFrom data.table := fread setnames
#' @import tidyr
#' @importFrom HDF5Array writeHDF5Array
#'
#'
### P-Value calculation
#'
#' @importFrom stats sd rbinom fisher.test na.omit
#' @importFrom VGAM rbetabinom vglm Coef pbetabinom pbetabinom.ab betabinomial
#'
### Miscelenious functions
#'
#' @importFrom BBmisc isScalarCharacter isScalarLogical chunk
#' @importFrom R.utils renameFile
#' @importFrom tools file_path_as_absolute
#' @importFrom methods as callNextMethod is new slot slot<- validObject
#' @importFrom utils browseURL capture.output
#'
#'

NULL

