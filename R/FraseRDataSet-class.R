#######################
## FraseRDataSet
## ====================

#' @include FraseRSettings-class.R
NULL

#' FraseRDataSet
#'
#' This class is designed to store the whole FraseR data set
#' needed for an analysis of a disease cohort
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
setClass("FraseRDataSet",
    slots = list(
        settings        = "FraseRSettings",
        splitReads      = "RangedSummarizedExperiment",
        nonSplicedReads = "RangedSummarizedExperiment"
    ),
    prototype = list(
        settings        = FraseRSettings(),
        splitReads      = SummarizedExperiment(rowRanges=GRanges()),
        nonSplicedReads = SummarizedExperiment(rowRanges=GRanges())
    )
)

## Validity
## ========

##
## Validating the correct type
##

.validateSettingsType <- function(object) {
    if(class(object@settings) != "FraseRSettings") {
        return("'settings' must be a FraseRSettings object.")
    }
    NULL
}

.validateSplitReadsType <- function(object) {
    if(class(object@splitReads) != "RangedSummarizedExperiment") {
        return("'splitReads' must be a RangedSummarizedExperiment object")
    }
    NULL
}

.validateNonSplicedReadsType <- function(object) {
    if(class(object@nonSplicedReads) != "RangedSummarizedExperiment") {
        return("'nonSplicedReads' must be a RangedSummarizedExperiment object")
    }
    NULL
}

## general validate function
.validateFraseRDataSet <- function(object) {
    c(
        .validateSettingsType(object),
        .validateSplitReadsType(object),
        .validateNonSplicedReadsType(object)
    )
}

setValidity2("FraseRDataSet", .validateFraseRDataSet)

## Constructor
## ==========

#' The constructor function for FraseRSettings
#'
#' @param ... Any parameters corresponding to the slots and their possible
#' values. See \linkS4class{FraseRDataSet}
#' @return A FraseRDataSet object.
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @examples
#'     fraser <- FraseRDataSet()
#'     fraser <- countRNAData(createTestFraseRSettings())
FraseRDataSet <- function(...) {
    return(new("FraseRDataSet", ...))
}

## Cosmetics
## =========

## show method for FraseRSettings
.showFraseRDataSet <- function(object) {

    # first show the setting part
    show(object@settings)

    cat("-------------------- Junction counts -------------------\n")
    show(object@splitReads)
    cat("\n\n")

    cat("-------------------- Splice site counts ----------------\n")
    show(object@nonSplicedReads)
    cat("\n\n")

}

setMethod("show", "FraseRDataSet", function(object) {
    .showFraseRDataSet(object)
})

## Getter and Setter
## =================

setGeneric("getDefaults", function(object, ...) standardGeneric("getDefaults"))

setMethod("getDefaults", "FraseRDataSet", function(object, what = NULL) {
    if(is.null(what))
        return(object)
    else
        return(slot(object, what))
})

setGeneric("setDefaults", function(object, ...) {
    standardGeneric("setDefaults")
})

setMethod("setDefaults", "FraseRDataSet", function(object, ...) {
    for(param in names(list(...))) {
        slot(object, param) <- list(...)[[param]]
    }
    return(object)
})



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
})

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
