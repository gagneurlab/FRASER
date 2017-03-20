#'
#'
#'
# save and load function for a fds dataset object
#setGeneric("saveFraseRDataSet",
#function(fds, dir=NULL, replace=FALSE, verbose=FALSE)
#    standardGeneric("saveFraseRDataSet"))
#setGeneric("loadFraseRDataSet",
#    function(dir) standardGeneric("loadFraseRDataSet"))
#'
#' #' @export
#' setMethod("saveFraseRDataSet", "FraseRDataSet",
#'     function (fds, dir=NULL, replace=FALSE, verbose=FALSE) {
#'         if(is.null(dir)){
#'             dir <- outputFolder(fds@settings)
#'         }
#'         outDir <- file.path(dir, "savedObjects")
#'         if(!dir.exists(outDir)){
#'             dir.create(outDir, recursive=TRUE)
#'         }
#'
#'         # convert assays to matrix to save them as HDF5Assay
#'         fds <- .assay2Matrix(fds)
#'
#'         # save spliceReads as HDF5
#'         fds@splitReads <- saveHDF5SummarizedExperiment(fds@splitReads,
#'                 dir=file.path(outDir, "splitReads"), replace=replace, verbose=verbose
#'         )
#'         # save nonSplicedReads as HDF5
#'         fds@nonSplicedReads <- saveHDF5SummarizedExperiment(
#'             fds@nonSplicedReads, dir=file.path(outDir, "nonSplicedReads"),
#'             replace=replace, verbose=verbose
#'         )
#'
#'         saveRDS(fds, file.path(outDir, "FraseRDataSet.RDS"))
#'         invisible(fds)
#'     }
#' )
#'
#' #' @export
#' setMethod("loadFraseRDataSet", "character", function (dir) {
#'   outDir <- file.path(dir, "savedObjects")
#'   if(!dir.exists(outDir)){
#'     stop(paste(
#'       "The given folder does not contain",
#'       "any saved FraseRDataSet objects."
#'     ))
#'   }
#'
#'   # load the main object
#'   fds <- readRDS(file.path(outDir, "FraseRDataSet.RDS"))
#'
#'   # load the SummarizedExperiment slots
#'   fds@splitReads <- loadHDF5SummarizedExperiment(
#'     dir=file.path(outDir, "splitReads")
#'   )
#'   fds@nonSplicedReads <- loadHDF5SummarizedExperiment(
#'     dir=file.path(outDir, "nonSplicedReads")
#'   )
#'
#'   return(fds)
#' })
#'
#' @noRd
.testinghdf5 <- function (dir = "my_h5_se"){
  outDir <- oldDir
  dir <- file.path(outDir, "splitReads")
  library(rhdf5)
  library(HDF5Array)
  if (!isSingleString(dir))
    stop(wmsg("'dir' must be a single string specifying the path ",
              "to the directory containing ", .THE_EXPECTED_STUFF))
  h5_path <- file.path(dir, "assays.h5")
  rds_path <- file.path(dir, "se.rds")
  if (!file.exists(h5_path) || !file.exists(rds_path))
    .stop_if_bad_dir(dir)
  h5_content <- try(rhdf5::h5ls(h5_path), silent = TRUE)
  if (inherits(h5_content, "try-error"))
    .stop_if_bad_dir(dir)
  h5_datasets <- h5_content[, "name"]
  ans <- readRDS(rds_path)
  if (!is(ans, "SummarizedExperiment"))
    .stop_if_bad_dir(dir)
  for (i in seq_along(assays(ans))) {
    a <- assay(ans, i, withDimnames = FALSE)
    # if (!is(a, "HDF5Array") || !identical(a@seed@file, "assays.h5") ||
    #     !(a@seed@name %in% h5_datasets))
    #   .stop_if_bad_dir(dir)
    a@seed@file <- file.path(dir, a@seed@file)
    assay(ans, i) <- a
  }
  str(assays(ans))
}


#'
#'
#' @export
saveFraseRDataSet <- function (fds, dir=NULL) {
    require(rhdf5)
    require(HDF5Array)

    if(is.null(dir)){
        dir <- outputFolder(fds@settings)
    }
    outDir <- file.path(dir, "savedObjects")
    if(!dir.exists(outDir)){
        dir.create(outDir, recursive=TRUE)
    }

    # over each se object
    for(type in c("splitReads", "nonSplicedReads")){
        assays <- assays(slot(fds, type))
        h5File <- .getFraseRHDF5File(fds, type, add=FALSE)
        # for each assay
        for(i in seq_along(assays)){
            assays <- .addAssayAsHDF5(assays, i, h5File)
        }
        assays(slot(fds, type)) <- assays
    }

    message("Writing final FraseR object.")
    saveRDS(fds, file.path(outDir, "fds-object.RDS"))

    if(endsWith(h5File, "1.h5")){
        for(type in c("splitReads", "nonSplicedReads")){
            unlink(file.path(outDir, paste0(type, ".h5")))
        }
    } else {
        for(type in c("splitReads", "nonSplicedReads")){
            f <- file.path(outDir, paste0(type, "1.h5"))
            if(file.exists(f)){
                unlink(f)
            }
        }
    }

    return(fds)
}

#'
#' creates the correct and needed HDF5 file
#'
#' @noRd
.getFraseRHDF5File <- function(fds, readType, add){
    dir <- outputFolder(fds@settings)
    outDir <- file.path(dir, "savedObjects")
    h5File <- file.path(outDir, paste0(readType, ".h5"))
    if(xor(file.exists(h5File), add)){
        h5File <- gsub(".h5", "1.h5", h5File)
    }
    return(h5File)
}

#'
#' save assay as hdf5 object and add it to the dataset
#'
#' @noRd
.addAssayAsHDF5 <- function(assays, idx, h5File){
    matrixGiven <- FALSE
    if(is.character(idx)){
        matrixGiven <- TRUE
        assays <- list(assays)
        aName  <- idx
        idx    <- 1
    } else {
        aName <- names(assays)[idx]
    }
    aMat <- as.matrix(as.data.table(assays[[idx]]))
    message("Writing assay: ", aName, " to file: ", h5File)
    ah5 <- writeHDF5Array(aMat, h5File, aName, verbose=TRUE)
    if(matrixGiven){
        return(ah5)
    }
    assays[[idx]] <- ah5
    return(assays)
}

#'
#'
#' @export
loadFraseRDataSet <- function(dir){
    require(rhdf5)
    require(HDF5Array)
    #dir <- file.path(outputFolder(fds@settings))
    #type <- "splitReads"
    #i <- 1

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
