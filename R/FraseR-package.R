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
#' @import BiocParallel
#'
#'
### GRange/Experiment packages
#'
#' @import GenomicAlignments
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom IRanges from
#' @importFrom IRanges to
#' @importFrom IRanges IRanges
#' @importFrom Rsubread featureCounts
#' @importFrom Rsamtools ScanBamParam
#'
#'
### Annotation
#'
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt getBM
#'
#'
### Plotting
#'
#' @importFrom graphics plot
#' @importFrom plotly plot_ly
#' @importFrom plotly subplot
#' @importFrom plotly layout
#' @importFrom htmlwidgets saveWidget
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom DelayedArray rowMeans
#' @importFrom rhdf5 h5ls
#'
#'
### Data handling
#'
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table fread
#' @import tidyr
#' @importFrom HDF5Array writeHDF5Array
#'
#'
### P-Value calculation
#'
#' @importFrom stats sd
#' @importFrom stats rbinom
#' @importFrom stats fisher.test
#' @importFrom stats na.omit
#' @importFrom VGAM rbetabinom
#' @importFrom VGAM vglm
#' @importFrom VGAM Coef
#' @importFrom VGAM pbetabinom
#' @importFrom VGAM betabinomial
#'
### Miscelenious functions
#'
#' @importFrom BBmisc isScalarCharacter
#' @importFrom BBmisc isScalarLogical
#' @importFrom R.utils renameFile
#' @importFrom tools file_path_as_absolute
#' @importFrom methods as
#' @importFrom methods callNextMethod
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom methods slot<-
#' @importFrom methods validObject
#' @importFrom utils browseURL
#' @importFrom utils capture.output
#'
#'

NULL

