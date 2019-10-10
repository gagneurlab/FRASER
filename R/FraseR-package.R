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
#' @importFrom plotly plot_ly subplot layout add_trace ggplotly
#' @importFrom gplots col2hex barplot2
#' @importFrom htmlwidgets saveWidget
#' @importFrom LSD heatscatter colorpalette
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
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
#' @importFrom rhdf5 h5ls H5Fopen H5Fclose
#'
### P-Value calculation
#'
#' @importFrom stats sd rbinom fisher.test na.omit p.adjust ppoints qbeta rnorm
#'          predict cor cutree dbinom dist hclust lm optim optimize pbinom
#'          plogis qlogis rlnorm rnbinom
#' @importFrom VGAM rbetabinom vglm Coef pbetabinom pbetabinom.ab betabinomial
#'          dbetabinom.ab dbetabinom
#'
### Miscelenious functions
#'
#' @importFrom BBmisc isScalarCharacter isScalarLogical chunk %nin%
#'          isScalarInteger isFALSE is.error isScalarValue isScalarNumeric
#' @importFrom R.utils renameFile withTimeout
#' @importFrom tools file_path_as_absolute
#' @importFrom methods as callNextMethod is new slot slot<- validObject
#' @importFrom utils browseURL capture.output sessionInfo
#'
#'
#
### To be added into the functions above
#
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom grDevices dev.off adjustcolor pdf colorRampPalette
#' @importFrom graphics abline axis grid legend lines title text points polygon
#'          hist
#' @importFrom plotly event_data
#' @importFrom GenomeInfoDb seqlevels<- seqlevels seqlengths
#'          keepStandardChromosomes
#' @importFrom shiny renderUI reactive renderPrint renderPlot
#' @importFrom DelayedArray rowMaxs rowMeans path<-
#' @importFrom data.table rbindlist
#' @importFrom DelayedMatrixStats rowMedians rowSds colMeans2 rowMeans2 rowQuantiles
#' @importFrom stats runif median quantile
#' @importFrom extraDistr rdirmnom dbbinom
#' @importFrom PRROC pr.curve
#' @importFrom ggplot2 ggtitle xlab ylab ggplot geom_point geom_line geom_smooth aes
#'          geom_line geom_hline geom_vline geom_abline geom_segment geom_ribbon
#'          scale_color_manual scale_x_log10 scale_y_log10 scale_color_gradientn
#'          labs theme_bw scale_color_discrete annotate geom_histogram theme
#'          scale_fill_manual
#' @importFrom ggpubr ggarrange
#'
#' @importFrom MASS kde2d bandwidth.nrd
#'
#' @importFrom keras custom_metric layer_input k_variable layer_lambda k_log
#'          layer_dense constraint_minmaxnorm regularizer_l2 get_weights
#'          set_weights keras_model k_exp k_mean optimizer_adam use_python
#'          callback_terminate_on_naan callback_early_stopping
#' @importFrom tensorflow install_tensorflow
#'
#' @importFrom tibble as_tibble
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


