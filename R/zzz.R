#'
#' The .onLoad function which is run during package loading
#' 
#' @noRd
#' 
.onLoad <- function(libname, pkgname){
    op <- options()
    op.fraser <- list(
        `FRASER-hdf5-chunk-nrow` = 30000,
        `FRASER-hdf5-chunk-ncol` = 20,
        `FRASER.pseudoCount` = 0.1,
        `FRASER.minSamplesForDelayed` = 1000,
        `FRASER.maxSamplesNoHDF5` = 20,
        `FRASER.maxJunctionsNoHDF5` = 1000)
    
    toset <- !(names(op.fraser) %in% names(op))
    if(any(toset)){
        options(op.fraser[toset])
    }
}
