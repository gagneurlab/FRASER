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
#' @import GenomicRanges
#' @import SummarizedExperiment
#  @import IRanges # not needed
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
#' @import plotly
#' 
#' 
### Data handling
#'
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table fread
#' @import tidyr
#' 
#' 
### P-Value calculation 
#' @importFrom VGAM rbetabinom
#' @importFrom VGAM vglm
#' @importFrom VGAM Coef
#' @importFrom VGAM pbetabinom
#' 
#' 
NULL
