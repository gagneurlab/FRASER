#'
#' The .onLoad function which is run during package loading
#' 
#' @noRd
#' 
.onLoad <- function(libname, pkgname){
    op <- options()
    op.fraser <- list(
        `FraseR-hdf5-chunk-nrow` = 30000,
        `FraseR-hdf5-chunk-ncol` = 20,
        `FraseR.pseudoCount` = 1,
        `FraseR.minSamplesForDelayed` = 1000)
    
    toset <- !(names(op.fraser) %in% names(op))
    if(any(toset)){
        options(op.fraser[toset])
    }
}
