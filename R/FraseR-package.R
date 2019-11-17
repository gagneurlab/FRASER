#'
#' FraseR: A package providing a workflow to detect aberrant splicing events
#'     in RNA sequencing data in a rare disease cohort
#'
### PACKAGE NAME
#'
#' @name FraseR
#'
#'
#' @import data.table
#'
#' @importFrom Biobase rowMax
#'
### Parallel computing
#'
#' @importFrom parallel mclapply detectCores
#' @import BiocParallel
#' @importFrom pcaMethods pca loadings
#'
### GRange/Experiment/bamFile packages
#' @importFrom BiocGenerics updateObject counts counts<- strand strand<-
#' @importFrom GenomicAlignments junctions readGAlignments summarizeJunctions
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames
#'          colData colData<- rowData rowRanges rowRanges<- SummarizedExperiment
#' @importFrom GenomicRanges findOverlaps granges GRanges GRangesList
#'          makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps from to IRanges ranges
#' @importFrom Rsamtools ScanBamParam scanBamHeader bamMapqFilter
#'          bamWhich bamWhich<-
#' @importFrom Rsubread featureCounts
#'
### Annotation
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#'
### Plotting
#'
#' @importFrom graphics abline axis box grid hist legend lines mtext par plot
#'          points polygon text title
#'
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
#' @importFrom tidyr %>%
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
#' @importFrom BBmisc isScalarCharacter isScalarLogical chunk %nin% seq_row
#'          seq_col isScalarInteger isFALSE is.error isScalarValue 
#'          isScalarNumeric
#' @importFrom R.utils renameFile withTimeout
#' @importFrom tools file_path_as_absolute
#' @importFrom methods as callNextMethod is new show slot slot<- validObject
#' @importFrom utils browseURL capture.output sessionInfo
#'
#'
#
### To be added into the functions above
#
#' @importFrom S4Vectors DataFrame metadata Rle SimpleList mcols mcols<-
#'          start end metadata metadata<-
#' @importFrom grDevices dev.off adjustcolor pdf colorRampPalette
#' @importFrom plotly event_data
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevels<- seqlevels
#'          seqlengths seqlengths<- seqlevelsStyle<- seqnames seqinfo
#' @importFrom shiny renderUI reactive renderPrint renderPlot
#' @importFrom DelayedArray rowMaxs rowMeans path<-
#' @importFrom DelayedMatrixStats colSds rowMedians rowSds colMeans2 rowMeans2
#'          rowQuantiles
#' @importFrom matrixStats colMaxs colMedians colMins
#' @importFrom stats runif median quantile
#' @importFrom extraDistr rdirmnom dbbinom pbbinom
#' @importFrom PRROC pr.curve
#' @importFrom ggplot2 ggtitle xlab ylab ggplot geom_point geom_line geom_smooth aes
#'          geom_line geom_hline geom_vline geom_abline geom_segment geom_ribbon
#'          scale_color_manual scale_x_log10 scale_y_log10 scale_color_gradientn
#'          labs theme_bw scale_color_brewer scale_color_discrete
#'          scale_linetype_manual annotate geom_histogram theme
#'          scale_fill_manual xlim scale_colour_manual
#' @importFrom ggpubr ggarrange
#'
#' @importFrom MASS kde2d bandwidth.nrd
#'
#' @importFrom keras callback_terminate_on_naan callback_early_stopping
#'          compile constraint_minmaxnorm custom_metric keras_model k_exp k_mean
#'          k_log k_set_value k_variable layer_input layer_lambda layer_dense
#'          regularizer_l2 get_weights set_weights optimizer_adam use_python
#'
#' @importFrom tensorflow install_tensorflow evaluate
#'
#' @importFrom tibble as_tibble
#'
#' @useDynLib FraseR
#'



NULL

#' TODO: Should be removed in the end!
#' @noRd
globalVariables(c(".", "J", ".N", ".asDataFrame", "End", "FN", "HTML", "Start", "TP",
        "deltaPsi", "curgr", "gene", "lty", "hgnc_symbol", "id",
        "ldat", "p.adj", "pval", "pvalue", "shinyFds", "shinyFdsRes",
        "sampleID", "sampleGroup", "chr", "symbol", "type", "pseudocount",
        "..sum_cols", "acceptorGroupID", "acceptorGroupSize", "Aggregation",
        "aroc", "donorGroupID", "donorGroupSize", "dpsi", "exMask", "fact",
        "fdrByFeature", "fds_inputK", "fds_inputN", "fds_new", "fds_rhoIn",
        "featureID", "features", "flog.debug", "flog.error", "flog.fatal",
        "flog.info", "flog.namespace", "flog.trace", "flog.warn", "geneID",
        "idx", "iteration", "Iteration", "IterSteps", "junctionID", "k", "Loss",
        "model", "mu", "n", ",nsubset", "o3", "o5", "obsPsi", "os", "pa",
        "padj", "passed", "pByFeature", "pointNr", "predPsi", "psi3", "psi5",
        "psiType", "psiValue", "seqlength", "seqlevel", "Step", "traceNr",
        "V1", "value", "zscore"),
        package="FraseR")
options("FraseR.pseudoCount"=1)
