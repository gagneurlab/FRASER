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
#' @param rewrite logical if the object should be rewritten. This makes sense if
#'             you have filtered or subsetted the object and want to save only
#'             the subsetted version
#'
#' @examples
#' fds <- countRNAData(createTestFraseRSettings())
#' fdsSaved <- saveFraseRDataSet(fds)
#' fdsSaved
#'
#' fdsLoaded <- loadFraseRDataSet(dir=workingDir(fds), name=name(fds))
#' fdsLoaded
#'
#' testthat::expect_equivalent(fdsSaved, fdsLoaded)
#'
#' @aliases loadFraseRDataSet saveFraseRDataSet
#' @rdname loadFraseRDataSet
#' @export
loadFraseRDataSet <- function(dir, name=NULL, upgrade=FALSE){
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
        } else {
            path(assay(fds, aname)) <- afile
        }
    }

    return(fds)
}


#' @rdname loadFraseRDataSet
#' @export
saveFraseRDataSet <- function(fds, dir=NULL, name=NULL, rewrite=FALSE) {

    # check input
    stopifnot(class(fds) == "FraseRDataSet")
    if(is.null(dir)) dir <- workingDir(fds)
    stopifnot(isScalarCharacter(dir))
    if(is.null(name)) name <- name(fds)
    stopifnot(isScalarCharacter(name))

    outDir <- file.path(dir, "savedObjects", nameNoSpace(name))
    if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

    # over each assay object
    for(aname in assayNames(fds)){
        assay <- assay(fds, aname)
        assay(fds, aname) <- saveAsHDF5(fds, aname, assay, rewrite=rewrite)
    }
    name(fds) <- name

    rdsFile <- file.path(outDir, "fds-object.RDS")
    message(date(), ": Writing final FraseR object ('", rdsFile, "').")
    saveRDS(fds, rdsFile)

    return(fds)
}


#'
#' saves the given assay as HDF5 array on disk
#' @noRd
saveAsHDF5 <- function(fds, name, object=NULL, rewrite=FALSE){
    if(is.null(object)) object <- assay(fds, name)

    if(isTRUE(dontWriteHDF5(fds))){
        message(date(), ": Dont save HDF5 for assay: ", name)
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
    message(date(), ": Preparing data for HDF5 conversion: ", name)
    aMat <- as(object, "matrix")
    message(date(), ": Writing data: ", name, " to file: ", h5File)
    h5 <- writeHDF5Array(aMat, h5FileTmp, name, verbose=FALSE, level=0)

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

