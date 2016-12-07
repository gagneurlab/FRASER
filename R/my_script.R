#' 
#' This is my first function and some doc
#'
#' @param arg this is a dummy parameter
#' @return NULL nothing
#' @export
#' @examples
#'   my_first_function()
#'
my_first_function <- function(arg = TRUE){
    message("working currently")
}

#'
#' Counting the split reads and the non spliced reads
#' in the given bamfile list. Can be computed in parallel
#'
#' @param bamfiles a list of bamfiles 
#' @param BPPARAM the BiocParallel param object to run jobs in parallel
#' @return SummarizedExperiment count data
#' @export
#' @examples
#'   count_rna("./data/my.bam")
#'
count_rna <- function(bamfiles, BPPARAM){
    message("counting rna files with BiocParallel")
}
