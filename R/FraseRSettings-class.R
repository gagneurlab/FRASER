#######################
## FraseRSettings
## ====================


#' FraseRSettings
#'
#' This class is designed to store the sample annotation 
#' and settings for the computation of the FraseR package
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
setClass("FraseRSettings",
        slots = list(
                sampleData = "data.table",
                bamParams = "ScanBamParam",
                parallel = "BiocParallelParam",
                method = "character",
                strandSpecific = "logical"
        ),
        prototype = list(
                sampleData = data.table(
                        sampleID="testSample", 
                        bamFile=system.file("extdata", "bam", "MUC1344.bam", package="FraseR")
                ),
                bamParams = Rsamtools::ScanBamParam(),
                parallel = BiocParallel::registered()[[1]],
                method = "zscore",
                strandSpecific = FALSE
        )
)

## Validity
## ========

##
## Validating the correct type
##

.validateSampleDataType <- function(object) {
    if(!any(class(object@sampleData) %in% "data.table")) {
        return("'sampleData' must be a data.frame or data.table object")
    }
    if(!any("sampleID" %in% colnames(object@sampleData))){
        return("The 'sampleData' table needs a sampleID column.")
    }
    if(any(object@sampleData[,duplicated(sampleID)])){
        return("The 'sampleID' column in the 'sampleData' needs to be unique.")
    }
    if(!any("bamFile" %in% colnames(object@sampleData))){
        return("The 'sampleData' table needs a bamFile column.")
    }
    NULL
}

.validateBamParamsType <- function(object) {
    if(class(object@bamParams) != "ScanBamParam") {
        return("'bamParams' must be a ScanBamParam object")
    }
    NULL
}

.validateMethodType <- function(object) {
    valid_methods <- c("zscore", "DESeq", "Martin")
    if(class(object@method) != "character" || !any(object@method %in% valid_methods)) {
        return(paste("'method' must be one of the following types: ", valid_methods))
    }
    NULL
}

.validateStrandSpecificType <- function(object) {
   if(class(object@strandSpecific) != "logical") {
        return(paste("'strandSpecific' must be of type logical."))
    }
    NULL
}

.validateParallelType <- function(object) {
    if(!is(object@parallel, "BiocParallelParam")) {
        return("'parallel' must be a BiocParallelParam object")
    }
    NULL
}

## general validate function
.validateFraseRSettings <- function(object) {
    c(
        .validateSampleDataType(object),
        .validateBamParamsType(object),
        .validateParallelType(object),
        .validateMethodType(object),
        .validateStrandSpecificType(object)
    )
}

setValidity2("FraseRSettings", .validateFraseRSettings)

## Constructor
## ==========

#' The constructor function for FraseRSettings
#' 
#' @param ... Any parameters corresponding to the slots and their possible
#' values. See \linkS4class{FraseRSettings}
#' @return A FraseRSettings object.
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @examples 
#'     settings <- FraseRSettings()
#'     settings <- FraseRSettings(fread(system.file(
#'             "extdata", "sampleTable.tsv", package="FraseR", mustWork=TRUE
#'     )))
FraseRSettings <- function(...) {
    return(new("FraseRSettings", ...)) 
}

## Cosmetics
## =========

## show method for FraseRSettings
.showFraseRSettings <- function(settings) {
    
    cat("-------------------- Sample data table -----------------\n")
    show(settings@sampleData)
    cat("\n\n")
    
    cat("-------------------- FraseR Settings --------------------\n")
    show(paste("Statistical method:", settings@method))
    show(paste("Strand specific data:", settings@strandSpecific))
    cat("\n\n")
    
    cat("-------------------- Parallel backend ------------------\n")
    show(settings@parallel)
    cat("\n\n")
    
    cat("-------------------- BAM parameters --------------------\n")
    show(settings@bamParams)
    cat("\n\n")
    
}

setMethod("show", "FraseRSettings", function(object) {
    .showFraseRSettings(object)
})

## Getter and Setter
## =================

setGeneric("getDefaults", function(object, ...) standardGeneric("getDefaults"))

setMethod("getDefaults", "FraseRSettings",
          function(object, what = NULL) {
              if(is.null(what)) return(object)
              else return(slot(object, what))
          })

setGeneric("setDefaults", function(object, ...) {
    standardGeneric("setDefaults")
})

setMethod("setDefaults", "FraseRSettings",
          function(object, ...) {
              for(param in names(list(...))) {
                  slot(object, param) <- list(...)[[param]]
                  if(param == "parallel") {
                      BiocParallel::register(slot(object, "parallel"))
                  }
              }
              return(object)
          })