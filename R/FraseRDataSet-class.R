#######################
## FraseRDataSet
## ====================


#' FraseRDataSet
#'
#' This class is designed to store the whole FraseR data set
#' needed for an analysis of a disease cohort
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
setClass("FraseRDataSet",
    contains="RangedSummarizedExperiment",
    slots = list(
        name            = "character",
        method          = "character",
        parallel        = "BiocParallelParam",
        bamParam        = "ScanBamParam",
        strandSpecific  = "logical",
        workingDir      = "character",
        nonSplicedReads = "RangedSummarizedExperiment"
    ),
    prototype = list(
        name            = "Data Analysis",
        method          = "betaBin",
        parallel        = SerialParam(),
        bamParam        = ScanBamParam(mapqFilter=10),
        strandSpecific  = FALSE,
        workingDir      = file.path(Sys.getenv("HOME"), "FraseR"),
        nonSplicedReads = SummarizedExperiment(rowRanges=GRanges())
    )
)

## Validity
## ========

#
# check sample annotation within the colData slot of the SE object
#
validateSampleAnnotation <- function(object) {
    sampleData <- as.data.table(colData(object))
    if(!"sampleID" %in% colnames(sampleData)){
        return("Please provide a 'sampleID' column with a ID for each sample.")
    }
    if(any(sampleData[,duplicated(sampleID)])){
        return("The 'sampleID' column needs to be unique.")
    }
    if(!any("bamFile" %in% colnames(sampleData))){
        return("Please provide a 'bamFile' column.")
    }
    NULL
}

validateName <- function(object){
    if(!isScalarCharacter(object@name)){
        return("The name of the experiment should be of type 'character'.")
    }
    if(length(object@name) == 0){
        return("The experiment name can not be empty.")
    }
    if(!grep("^[a-zA-Z0-9 ._-]+$", object@name, perl=TRUE)){
        return(paste("For readabilty the name of the experiment should only ",
                "contain the following characters: 'a-zA-Z0-9 ._-'"
        ))
    }
    NULL
}

validateMethod <- function(object) {
    validMethods <- c("Fisher", "betaBin", "DESeq2", "Martin")
    if(!isScalarCharacter(object@method) || !object@method %in% validMethods) {
        return(paste0("The selected method must be one of the following: ",
                paste(validMethods, collapse=", "), "."
        ))
    }
    NULL
}

validateParallel <- function(object) {
    if(!is(object@parallel, "BiocParallelParam")) {
        return("The 'parallel' option must be a BiocParallelParam object.")
    }
    NULL
}

validateBamParam <- function(object) {
    if(class(scanBamParam(object)) != "ScanBamParam") {
        return("The 'bamParam' option must be a ScanBamParam object.")
    }
    NULL
}

validateStrandSpecific <- function(object) {
    if(!isScalarLogical(object@strandSpecific)) {
        return(paste("The 'strandSpecific' option must be TRUE or FALSE."))
    }
    NULL
}


validateWorkingDir <- function(object) {
    if(!isScalarCharacter(object@workingDir)){
        return(paste("The path to the working directory needs",
                "to be a set as a character."
        ))
    }
    if(object@workingDir == ""){
        return("The working directory can not be empty.")
    }
    if(!dir.exists(object@workingDir)){
        message(date(), ": The given working directory '", object@workingDir,
                "' does not exists. We will create it."
        )
        dir.create(object@workingDir, recursive = TRUE)
    }
    if(file.access(object@workingDir, mode = 2) != 0){
        return(paste("Make sure we can write to the given working directory '",
                object@workingDir, "'."
        ))
    }
    NULL
}

validateNonSplicedReadsType <- function(object) {
    if(class(object@nonSplicedReads) != "RangedSummarizedExperiment") {
        return("'nonSplicedReads' must be a RangedSummarizedExperiment object")
    }
    if(length(object) != 0 && dim(object@nonSplicedReads)[2] != dim(object)[2]){
        return("The NSR dimensions are not correct. This is a internal error!")
    }
    ans <- validObject(object@nonSplicedReads)
    if(!isScalarLogical(ans) || ans == FALSE){
        return(ans)
    }
    NULL
}


## general validate function
validateFraseRDataSet <- function(object) {
    c(
        validateSampleAnnotation(object),
        validateName(object),
        validateMethod(object),
        validateParallel(object),
        validateBamParam(object),
        validateStrandSpecific(object),
        validateWorkingDir(object),
        validateNonSplicedReadsType(object)
    )
}
setValidity("FraseRDataSet", validateFraseRDataSet)


## Cosmetics (the show function)
## =============================

## show method for FraseRDataSet
showFraseRDataSet <- function(object) {
    # taken from SummarizedExperiment show function
    scat <- function(fmt, vals=character(), exdent=2, ...){
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("-------------------- Sample data table -----------------\n")
    sampleData <- as.data.table(colData(object))
    if(all(sapply(sampleData[,bamFile], isScalarCharacter))){
        sampleData[,bamFile:=sapply(bamFile, function(str){
                if(nchar(str) <= 29) return(str)
                paste("...", substr(str, nchar(str) - 25, nchar(str)))
        })]
        sampleData[,bamFile:=gsub("... [^/]+/", ".../", bamFile)]
    }
    show(sampleData)
    cat("\n\n")
    if(length(object) > 0){
        cat(paste0("Number of samples:      ", dim(object)[2]), "\n")
        cat(paste0("Number of junctions:    ", length(object)), "\n")
        cat(paste0("Number of splice sites: ", length(nonSplicedReads(object))), "\n")
        scat("assays(%d):    %s\n", assayNames(object))
        cat("\n\n")
    }

    cat("----------------------- Settings -----------------------\n")
    cat(paste0("Analysis name:               ", name(object)), "\n")
    cat(paste0("Statistical method:          ", method(object)), "\n")
    cat(paste0("Analysis is strand specific: ", strandSpecific(object)), "\n")
    cat(paste0("Working directory:           '", workingDir(object), "'"), "\n")
    cat("\n\n")

    cat("-------------------- Parallel backend ------------------\n")
    show(parallel(object))
    cat("\n\n")

    cat("-------------------- BAM parameters --------------------\n")
    show(scanBamParam(object))
    cat("\n\n")
}

setMethod("show", "FraseRDataSet", function(object) {
    showFraseRDataSet(object)
})

## Constructor
## ==========

#'
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
FraseRDataSet <- function(colData=NULL, ...) {
    if(!is.null(colData)){
        if(is.data.table(colData)){
            colData <- DataFrame(colData)
        }
        return(new("FraseRDataSet", colData=colData, ...))
    }
    return(new("FraseRDataSet", ...))
}



