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
        parallel        = "BiocParallelParam",
        bamParam        = "ScanBamParam",
        strandSpecific  = "logical",
        workingDir      = "character",
        nonSplicedReads = "RangedSummarizedExperiment"
    ),
    prototype = list(
        name            = "Data Analysis",
        parallel        = SerialParam(),
        bamParam        = ScanBamParam(mapqFilter=0),
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
    if(any(duplicated(sampleData$sampleID))){
        return("The 'sampleID' column needs to be unique.")
    }
    if(!any("bamFile" %in% colnames(sampleData))){
        return("Please provide a 'bamFile' column.")
    }
    if(any(samples(object) != rownames(colData(object)))){
        return("Please set the rownames of your colData to the sampleIDs")
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
        return(paste("The 'strandSpecific' option must be 0L, 1L or 2L."))
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
    if(file.access(object@workingDir, mode = 4) != 0){
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

validateAssays <- function(object){
    if(length(assayNames(object)) > 1){
        if(any(duplicated(assayNames(object)))){
            return(
                "FraseR enforces unique assay names! Please provie such names."
            )
        }
    }
    NULL
}


## general validate function
validateFraseRDataSet <- function(object) {
    c(
        validateSampleAnnotation(object),
        validateName(object),
        validateParallel(object),
        validateBamParam(object),
        validateStrandSpecific(object),
        validateWorkingDir(object),
        validateNonSplicedReadsType(object),
        validateAssays(object)
    )
}
setValidity("FraseRDataSet", validateFraseRDataSet)


## Cosmetics (the show function)
## =============================

## show method for FraseRDataSet
showFraseRDataSet <- function(object) {
    if(dim(object)[2] < 1){
        cat("This Fraser object does not contain any sample! Please add one.")
        return()
    }

    # taken from SummarizedExperiment show function
    scat <- function(fmt, vals=character(), exdent=2, ...){
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("-------------------- Sample data table -----------------\n")
    sampleData <- as.data.table(colData(object))
    if(all(sapply(sampleData$bamFile, isScalarCharacter))){
        sampleData$bamFile <- gsub("... [^/]+/", ".../",
            sapply(sampleData$bamFile, function(str){
                if(nchar(str) <= 29) return(str)
                paste("...", substr(str, nchar(str) - 25, nchar(str)))
            })
        )
    }
    show(as_tibble(sampleData))
    cat("\n")
    if(length(object) > 0){
        cat(paste0("Number of samples:      ", dim(object)[2]), "\n")
        cat(paste0("Number of junctions:    ", length(object)), "\n")
        cat(paste0("Number of splice sites: ", length(nonSplicedReads(object))), "\n")
        scat("assays(%d):    %s\n", assayNames(object))
        cat("\n")
    }

    cat("----------------------- Settings -----------------------\n")
    cat(paste0("Analysis name:               ", name(object)), "\n")
    cat(paste0("Analysis is strand specific: ", strandSpecific(object)), "\n")
    cat(paste0("Working directory:           '", workingDir(object), "'"), "\n")
    cat("\n")

    cat("-------------------- Parallel backend ------------------\n")
    # show(parallel(object))
    cat(paste0("Type: ", as.character(class(parallel(object))),
            "\tWorkers: ", bpworkers(parallel(object)),
            "\tTasks: ", bptasks(parallel(object))
    ))
    cat("\n\n")

    cat("-------------------- BAM parameters --------------------\n")
    if(identical(scanBamParam(FraseRDataSet()), scanBamParam(object))){
        cat(paste0("Default used with: ",
                "bamMapqFilter=", bamMapqFilter(scanBamParam(object))
        ))
    } else {
        show(scanBamParam(object))
    }
    cat("\n\n")
}

setMethod("show", "FraseRDataSet", function(object) {
    showFraseRDataSet(object)
})

## Constructor
## ==========

#'
#' The FraseR dataset object
#'
#' Constructs an FraseR object based on the given input. It can take only the
#' annotation (colData) or count tables (junctions/spliceSites).
#'
#' @param colData A DataFrame containing the annotation of the samples
#' @param junctions A matrix like object containing the raw counts for each
#'                junction. It requires the \code{start} and the \code{endID}
#'                column that identifies the corresponding splice site for
#'                the given junction.
#' @param spliceSites A matrix like object containing the raw counts for each
#'                splice site. it requires the \code{spliceSiteID} and the
#'                \code{type} column that gives the ID and the type of the
#'                given splice site. The ID maps back to the junction.
#' @param ... Any parameters corresponding to the slots and their possible
#'                values. See \linkS4class{FraseRDataSet}
#' @return A FraseRDataSet object.
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @examples
#'     fraser <- FraseRDataSet()
#'     fraser <- countRNAData(createTestFraseRSettings())
FraseRDataSet <- function(colData=NULL, junctions=NULL, spliceSites=NULL, ...) {
    if(!is.null(colData)){
        if(is.data.table(colData)){
            colData <- DataFrame(colData)
        }
        if(is.null(rownames(colData))){
            rownames(colData) <- colData[['sampleID']]
        }
        if(is.null(junctions) & is.null(spliceSites))
            return(new("FraseRDataSet", colData=colData, ...))
        if(is.null(junctions)){
            stop("Please provide junctions counts if you provide ",
                    "spliceSite counts.")
        }
        if(is.null(spliceSites)){
            stop("Please provdie splice site counts if you provide ",
                    "junction counts.")
        }
        if(is.data.frame(junctions)){
            junctions <- makeGRangesFromDataFrame(junctions,
                    keep.extra.columns=TRUE)
        }
        if(is.data.frame(spliceSites)){
            spliceSites <- makeGRangesFromDataFrame(spliceSites,
                    keep.extra.columns=TRUE)
        }
        nsr <- SummarizedExperiment(
                rowRanges=spliceSites[,c("spliceSiteID", "type")],
                assays=SimpleList(rawCountsSS=as.data.frame(
                        mcols(spliceSites)[colData[,"sampleID"]]), a=NULL)[1])
        se <- SummarizedExperiment(
                rowRanges=junctions[,c("startID", "endID")],
                colData=colData,
                assays=SimpleList(rawCountsJ=as.data.frame(
                        mcols(junctions)[colData[,"sampleID"]]), a=NULL)[1])
        return(new("FraseRDataSet", se, nonSplicedReads=nsr, ...))
    }
    return(new("FraseRDataSet", ...))
}



