#'
#' Load a saved FraseR object into memory
#'
#' @param dir a path to the working directory of FraseR
#'
#' @example
#'   loadFraseRDataSet('~/FraseR')
#'
#' @export
loadFraseRDataSet <- function(dir){
    if(is.null(dir)){
        stop("dir: canot be null")
    }
    if(!dir.exists(dir)){
        stop("The given dir does not exists: ", dir)
    }
    outDir <- file.path(dir, "savedObjects")

    fds <- readRDS(file.path(outDir, "fds-object.RDS"))
    return(fds)
}


#'
#' Saves the FraseRDataSet object on disk under the given working dir.
#' It furthermore uses HDF5 to save all internal assays
#'
#' After saving this can be loaded again with the corresnonding load function
#'
#' @param fds A FraseRDataSet object ot be saved
#' @param dir a directory name where to save the objects
#'             (replaces the working directory)
#' @examples
#'     fds <- countRNAData(createTestFraseRSettings())
#'     saveFraseRDataSet(fds)
#'
#' @export
saveFraseRDataSet <- function(fds, dir=NULL) {

    # check input
    stopifnot(class(fds) == "FraseRDataSet")
    if(any(assayNames(fds) == "")) stop("Name of an assay can not be empty!")
    if(any(duplicated(assayNames(fds)))) stop("Assay names need to be unique!")

    if(is.null(dir)){
        dir <- workingDir(fds)
    }
    outDir <- file.path(dir, "savedObjects")
    if(!dir.exists(outDir)){
        dir.create(outDir, recursive=TRUE)
    }

    # over each se object
    assays <- assays(fds)
    for(aname in names(assays)){
        assay <- assay(fds, aname)
        assays[[aname]] <- saveAsHDF5(fds, aname, assay)
    }
    assays(fds) <- assays

    message(date(), ": Writing final FraseR object.")
    saveRDS(fds, file.path(outDir, "fds-object.RDS"))

    return(fds)
}


#'
#' saves the given assay as HDF5 array on disk
#' @noRd
saveAsHDF5 <- function(fds, name, object=NULL){
    if(is.null(object)){
        object <- assay(fds, name)
    }
    h5File <- getFraseRHDF5File(fds, name, add=FALSE)
    h5FileTmp <- gsub("$", ".save.tmp", h5File)
    if(file.exists(h5FileTmp)) unlink(h5FileTmp)

    # write new HDF5 data
    message(date(), ": Writing data: ", name, " to file: ", h5File)
    aMat <- as.matrix(as.data.table(object))
    h5 <- writeHDF5Array(aMat, h5FileTmp, name, verbose=TRUE)

    # override old h5 file if present and move tmp to correct place
    if(file.exists(h5File)) unlink(h5File)
    renameFile(h5FileTmp, h5File)
    h5@seed@file <- h5File

    return(h5)
}


#'
#' creates the correct and needed HDF5 file
#' @noRd
getFraseRHDF5File <- function(fds, aname, add){
    dir <- workingDir(fds)
    outDir <- file.path(dir, "savedObjects")
    h5File <- file.path(outDir, paste0(aname, ".h5"))
    return(h5File)
}

