#'
#' Load a saved FraseR object into memory
#'
#' @param dir a path to the working directory of FraseR
#'
#' @examples
#'   loadFraseRDataSet(file.path(Sys.getenv("HOME"), "FraseR"))
#'
#' @export
loadFraseRDataSet <- function(dir, name=NULL){
    # check dir
    if(is.null(dir)) stop("dir: can not be NULL")
    if(!isScalarCharacter(dir)) stop("dir: needs to be a character path name.")
    if(!dir.exists(dir)) stop("The given dir does not exists: ", dir)

    # check name
    if(is.null(name)) name <- "Data Analysis"
    if(!isScalarCharacter(name)) stop("name: needs to be a character dir name.")
    outDir <-file.path(dir, "savedObjects", nameNoSpace(name))
    if(!dir.exists(outDir)) stop("The analysis name does not exists: ", name)

    fds <- readRDS(gsub("//+", "/", file.path(outDir, "fds-object.RDS")))

    # set working dir and name correct
    workingDir(fds) <- dir
    name(fds) <- name

    # set the correct path of the assay seed file (if folder changed)
    for(obj in c(fds, nonSplicedReads(fds))){
        for(aname in names(obj@assays$data@listData)){
            afile <- getFraseRHDF5File(fds, aname)
            obj@assays$data@listData[[aname]]@seed@file <- afile
        }
    }

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
    if(is.null(dir)) dir <- workingDir(fds)
    stopifnot(isScalarCharacter(dir))

    outDir <- file.path(dir, "savedObjects", nameNoSpace(fds))
    if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # over each assay object
    assays <- assays(fds)
    for(aname in names(assays)){
        assay <- assay(fds, aname)
        assays[[aname]] <- saveAsHDF5(fds, aname, assay)
    }
    assays(fds) <- assays

    rdsFile <- file.path(outDir, "fds-object.RDS")
    message(date(), ": Writing final FraseR object ('", rdsFile, "').")
    saveRDS(fds, rdsFile)

    return(fds)
}


#'
#' saves the given assay as HDF5 array on disk
#' @noRd
saveAsHDF5 <- function(fds, name, object=NULL){
    if(is.null(object)) object <- assay(fds, name)

    h5File <- getFraseRHDF5File(fds, name)
    h5FileTmp <- gsub("$", ".save.tmp", h5File)
    if(file.exists(h5FileTmp)) unlink(h5FileTmp)

    # dont rewrite it if already there
    if("DelayedMatrix" %in% is(object) && object@seed@file == h5File){
        return(object)
    }

    # write new HDF5 data
    message(date(), ": Preparing data for HDF5 conversion: ", name)
    aMat <- as.matrix(as.data.table(object))
    message(date(), ": Writing data: ", name, " to file: ", h5File)
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
getFraseRHDF5File <- function(fds, aname){
    dir <- workingDir(fds)
    outDir <- file.path(dir, "savedObjects", nameNoSpace(fds))
    if(!dir.exists(outDir)) {
        dir.create(outDir, recursive=TRUE)
    }
    h5File <- file.path(file_path_as_absolute(outDir), paste0(aname, ".h5"))
    return(h5File)
}

