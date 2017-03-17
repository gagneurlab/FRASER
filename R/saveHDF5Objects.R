
# save and load function for a fds dataset object
setGeneric("saveFraseRDataSet",
function(fds, dir=NULL, replace=FALSE, verbose=FALSE)
    standardGeneric("saveFraseRDataSet"))
setGeneric("loadFraseRDataSet",
    function(dir) standardGeneric("loadFraseRDataSet"))

#' @export
setMethod("saveFraseRDataSet", "FraseRDataSet",
    function (fds, dir=NULL, replace=FALSE, verbose=FALSE) {
        if(is.null(dir)){
            dir <- outputFolder(fds@settings)
        }
        outDir <- file.path(dir, "savedObjects")
        if(!dir.exists(outDir)){
            dir.create(outDir, recursive=TRUE)
        }

        # convert assays to matrix to save them as HDF5Assay
        fds <- .assay2Matrix(fds)

        # save spliceReads as HDF5
        fds@splitReads <- saveHDF5SummarizedExperiment(fds@splitReads,
                dir=file.path(outDir, "splitReads"), replace=replace, verbose=verbose
        )
        # save nonSplicedReads as HDF5
        fds@nonSplicedReads <- saveHDF5SummarizedExperiment(
            fds@nonSplicedReads, dir=file.path(outDir, "nonSplicedReads"),
            replace=replace, verbose=verbose
        )

        saveRDS(fds, file.path(outDir, "FraseRDataSet.RDS"))
        invisible(fds)
    }
)

#' @export
setMethod("loadFraseRDataSet", "character", function (dir) {
  outDir <- file.path(dir, "savedObjects")
  if(!dir.exists(outDir)){
    stop(paste(
      "The given folder does not contain",
      "any saved FraseRDataSet objects."
    ))
  }

  # load the main object
  fds <- readRDS(file.path(outDir, "FraseRDataSet.RDS"))

  # load the SummarizedExperiment slots
  fds@splitReads <- loadHDF5SummarizedExperiment(
    dir=file.path(outDir, "splitReads")
  )
  fds@nonSplicedReads <- loadHDF5SummarizedExperiment(
    dir=file.path(outDir, "nonSplicedReads")
  )

  return(fds)
})

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
