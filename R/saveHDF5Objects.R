#'
#' Loading/Saving FraseRDataSets
#'
#' This is a convinient function to load and save a FraseRDataSet object.
#' It looks and saves the FraseRDataSet objects and HDF5 files on disk under
#' the given working dir. Internally it uses HDF5 files for all assays.
#'
#' @param fds A FraseRDataSet object ot be saved
#' @param dir A path where to save the objects (replaces the working directory)
#' @param name The analysis name of the project (saved within the `dir`)
#' @param file The file path to the fds-object.RDS file that should be loaded.
#' @param rewrite logical if the object should be rewritten. This makes sense if
#'             you have filtered or subsetted the object and want to save only
#'             the subsetted version
#' @param upgrade Should the version of the loaded object be updated?
#'
#' @examples
#' fds <- createTestFraseRSettings()
#' name(fds) <- "saveing_test"
#' 
#' # make sure the object is saved to disc
#' dontWriteHDF5(fds) <- FALSE
#' fdsSaved <- saveFraseRDataSet(fds)
#' fdsSaved
#' 
#' # load object from disc
#' fdsLoaded <- loadFraseRDataSet(dir=workingDir(fds), name=name(fds))
#' fdsLoaded
#'
#' all.equal(fdsSaved, fdsLoaded)
#' 
#' @return FraseRDataSet
#' @aliases loadFraseRDataSet saveFraseRDataSet
#' @rdname loadFraseRDataSet
#' @export
loadFraseRDataSet <- function(dir, name=NULL, file=NULL, upgrade=FALSE){
    # check if file is provided
    if(!is.null(file)){
        if(!missing(dir) | !is.null(name)){
            stop("You can only provide 'file' or 'dir' + 'name', ",
                    "but not all together.")
        }
        if(!file.exists(file) | !grepl("\\.(RDS|h5)$", file, perl=TRUE)){
            stop("Please provide the `fds-object.RDS` file.")
        }
        name <- basename(dirname(file))
        dir  <- dirname(dirname(dirname(file)))
    }
    
    # check dir
    if(is.null(dir)) stop("dir: can not be NULL")
    if(!isScalarCharacter(dir)) stop("dir: needs to be a character path name.")
    if(!dir.exists(dir)) stop("The given dir does not exists: ", dir)
    
    # check name
    if(is.null(name)) name <- "Data Analysis"
    if(!isScalarCharacter(name)) stop("name: needs to be a character dir name.")
    outDir <-file.path(dir, "savedObjects", nameNoSpace(name))
    if(!dir.exists(outDir)) stop("The analysis name does not exists: ", name)

    fdsFile <- gsub("//+", "/", file.path(outDir, "fds-object.RDS"))
    if(!file.exists(fdsFile)) {
        stop(paste("The main RDS file", fdsFile,
                "does not exists. Did you saved it correctly before?"
        ))
    }
    fds <- readRDS(fdsFile)
    
    # needs to be here due to our FraseR -> FRASER package change.
    # can be removed later if the full pipeline is rerun
    attributes(fds)$class <- structure("FraseRDataSet", package="FRASER")
    
    e <- try(assays(fds), silent=TRUE)
    if(is.error(e)){
        if(grepl("DelayedMatrix .* representation .* Please update it ",
                as.character(e))){
            if(isTRUE(upgrade)){
                for(i in assayNames(a)){
                    obj <- updateObject(a[[i]], verbose=TRUE)
                    a$data[[i]] <- obj
                }
            } else {
                stop(paste('Please upgrade the DelayedMatrix',
                        'objects or set the update option to TRUE\n\n'), e)
            }
        }
    }
    
    # set working dir and name correct
    workingDir(fds) <- dir
    name(fds) <- name
    
    # set the correct path of the assay seed file (if folder changed)
    for(aname in assayNames(fds)){
        message("Loading assay: ", aname)
        afile <- getFraseRHDF5File(fds, aname)
        if(!file.exists(afile)){
            warning(paste("Can not find assay file: ", aname, ".",
                    "The assay will be removed from the object."))
            assay(fds, aname) <- NULL
        } else if(afile != path(assay(fds, aname))) {
            path(assay(fds, aname)) <- afile
        }
    }

    return(fds)
}


#' @rdname loadFraseRDataSet
#' @export
saveFraseRDataSet <- function(fds, dir=NULL, name=NULL, rewrite=FALSE) {

    if(isTRUE(dontWriteHDF5(fds))){
        if(verbose(fds) > 0){
            message(date(), ": Dont save fds object.")
        }
        return(fds)
    }

    # check input
    stopifnot(is(fds, "FraseRDataSet"))
    if(is.null(dir)) dir <- workingDir(fds)
    stopifnot(isScalarCharacter(dir))
    if(is.null(name)) name <- name(fds)
    stopifnot(isScalarCharacter(name))

    outDir <- file.path(dir, "savedObjects", nameNoSpace(name))
    if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # over each assay object
    name(fds) <- name
    for(aname in assayNames(fds)){
        assay <- assay(fds, aname)
        assay(fds, aname) <- saveAsHDF5(fds, aname, assay, rewrite=rewrite)
    }

    rdsFile <- file.path(outDir, "fds-object.RDS")
    message(date(), ": Writing final FraseR object ('", rdsFile, "').")
    saveRDS(fds, rdsFile)

    return(fds)
}

#'
#' Defaults for HDF5 global options
#' To save the underlying assays
#' @noRd
options("FraseR-hdf5-chunk-nrow"=30000)
options("FraseR-hdf5-chunk-ncol"=20)

#'
#' saves the given assay as HDF5 array on disk
#' @noRd
saveAsHDF5 <- function(fds, name, object=NULL, rewrite=FALSE){
    if(is.null(object)) object <- assay(fds, name)

    # get defind chunk sizes
    chunkDims <- c(
        min(nrow(object), options()[['FraseR-hdf5-chunk-nrow']]),
        min(ncol(object), options()[['FraseR-hdf5-chunk-ncol']]))
    
    if(isTRUE(dontWriteHDF5(fds))){
        if(verbose(fds) > 3){
            message(date(), ": Dont save HDF5 for assay: ", name)
        }
        return(object)
    }
    h5File <- getFraseRHDF5File(fds, name)
    h5FileTmp <- paste0(h5File, ".", as.integer(abs(rnorm(1))*100), ".save.tmp")
    if(file.exists(h5FileTmp)) unlink(h5FileTmp)

    # dont rewrite it if already there
    if(!rewrite && "DelayedMatrix" %in% is(object) && path(object) == h5File){
        return(object)
    }

    # write new HDF5 data
    if(verbose(fds) > 2) {
        message(date(), ": Preparing data for HDF5 conversion: ", name)
    }
    aMat <- as.matrix(object)
    #aMat <- as(object, "matrix") # this command breaks the code in some cases
    if(verbose(fds) > 1) {
        message(date(), ": Writing data: ", name, " to file: ", h5File)
    }
    h5 <- writeHDF5Array(aMat, h5FileTmp, name, verbose=FALSE, 
            chunkdim=chunkDims)

    # override old h5 file if present and move tmp to correct place
    if(file.exists(h5File)) unlink(h5File)
    renameFile(h5FileTmp, h5File)
    path(h5) <- h5File

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

