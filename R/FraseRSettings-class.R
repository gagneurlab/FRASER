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
        strandSpecific = "logical",
        outputFolder = "character"
    ),
    prototype = list(
        sampleData = data.table(
            sampleID="testSample", 
            bamFile=system.file(
                "extdata", "bam", "sample1.bam", package="FraseR"
            )
        ),
        bamParams = Rsamtools::ScanBamParam(),
        parallel = BiocParallel::registered()[[1]],
        method = "betaBin",
        strandSpecific = FALSE,
        outputFolder = file.path(Sys.getenv("HOME"), "FraseR")
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
    valid_methods <- c("Fisher", "betaBin", "DESeq2", "Martin")
    if(class(object@method) != "character" || 
        !any(object@method %in% valid_methods)) {
        return(paste(
            "'method' must be one of the following types: ", valid_methods
        ))
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

.validateOutputFolder <- function(object) {
    if(is.null(object@outputFolder)){
        return(NULL)
    }
    if(class(object@outputFolder) != "character"){
        return("'outputFolder' must be a valid path or NULL.")
    }
    if(object@outputFolder == ""){
        return("'outputFolder' should not be a empty string!")
    }
    if(!dir.exists(object@outputFolder)){
        dir.create(object@outputFolder, recursive = TRUE)
    }
    if(file.access(object@outputFolder, mode = 2) != 0){
        return("Make sure we can write to the 'outputFolder' or set it to NULL")
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
        .validateStrandSpecificType(object),
        .validateOutputFolder(object)
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
    show(paste0("Output folder for plots and cache: '", 
            settings@outputFolder, "'")
    )
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
    }
)

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
    }
)

# sampleData functions
setGeneric("sampleData", function(object) standardGeneric("sampleData"))
setGeneric("sampleData<-", signature = "object", 
           function(object, value) standardGeneric("sampleData<-"))

# sampleGroup functions
setGeneric("sampleGroup", function(object) standardGeneric("sampleGroup"))
setGeneric("sampleGroup<-", signature = "object", 
    function(object, value) standardGeneric("sampleGroup<-"))

# outputFolder functions
setGeneric("outputFolder", function(object) standardGeneric("outputFolder"))
setGeneric("outputFolder<-", signature = "object", 
   function(object, value) standardGeneric("outputFolder<-"))

# scanBamParam functions
setGeneric("scanBamParam", function(object) standardGeneric("scanBamParam"))
setGeneric("scanBamParam<-", signature = "object", 
           function(object, value) standardGeneric("scanBamParam<-"))

# parallel functions
setGeneric("parallel", function(object) standardGeneric("parallel"))
setGeneric("parallel<-", signature = "object", 
           function(object, value) standardGeneric("parallel<-"))

# method functions
setGeneric("method", function(object) standardGeneric("method"))
setGeneric("method<-", signature = "object", 
           function(object, value) standardGeneric("method<-"))

# samples functions
setGeneric("samples", function(object) standardGeneric("samples"))
setGeneric("samples<-", signature = "object", 
           function(object, value) standardGeneric("samples<-"))

# samples functions
setGeneric("bamFiles", function(object) standardGeneric("bamFiles"))
setGeneric("bamFiles<-", signature = "object", 
           function(object, value) standardGeneric("bamFiles<-"))


#'
#'  Get/Set the sampleData table the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A \code{data.table} with all sample infromation (id, bamfile, group)
#' @examples
#' settings <- createTestFraseRSettings()
#' sampleData(settings)
#' sampleData(settings) <- 1:length(settings)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname sampleData
setMethod("sampleData", "FraseRSettings", function(object) {
    return(slot(object, "sampleData"))
})

#' @export
#' @rdname sampleData
setReplaceMethod("sampleData", "FraseRSettings", function(object, value) {
    slot(object, "sampleData") <- value
    return(object)
})


#'
#'  Get/Set the sampleIDs from the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A \code{vector} with all sample IDs
#' @examples
#' settings <- createTestFraseRSettings()
#' samples(settings)
#' samples(settings) <- 1:length(settings)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname samples
setMethod("samples", "FraseRSettings", function(object) {
    return(sampleData(object)[,sampleID])
})

#' @export
#' @rdname samples
setReplaceMethod("samples", "FraseRSettings", function(object, value) {
    sampleData(object)[,sampleID:=value]
    return(object)
})


#'
#'  Get the group definition per sample based on the 
#' \code{sampleData} table slot
#' 
#' @param object A FraseRSettings object.
#' @return A \code{vector} with the group identifiers per sample object
#' @examples
#' settings <- createTestFraseRSettings()
#' sampleGroup(settings)
#' sampleGroup(settings) <- 1:length(settings)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname sampleGroup
setMethod("sampleGroup", "FraseRSettings", function(object) {
    data <- sampleData(object)
    if("group" %in% colnames(data)){
        return(data[,group])
    } else {
        return(data[,sampleID])
    }
})

#' @export
#' @rdname sampleGroup
setReplaceMethod("sampleGroup", "FraseRSettings", function(object, value) {
    sampleData(object)[,group:=value]
    return(object)
})


#'
#'  Get/Set the bamFiles from the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A \code{vector} with the bamFiles for each sample
#' @examples
#' settings <- createTestFraseRSettings()
#' bamFiles(settings)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname bamFiles
setMethod("bamFiles", "FraseRSettings", function(object) {
    bamFiles <- sampleData(object)[,bamFile]
    if(all(sapply(bamFiles, class) == "character")){
        bamFiles <- unlist(BamFileList(bamFiles))
    }
    return(bamFiles)
})

#' @export
#' @rdname bamFiles
setReplaceMethod("bamFiles", "FraseRSettings", function(object, value) {
    sampleData(object)[,bamFile:=value]
    return(object)
})


#'
#'  Get/Set the parallel param object from the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A parallel param object
#' @examples
#' settings <- createTestFraseRSettings()
#' parallel(settings)
#' parallel(settings) <- SerialParam()
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname parallel
setMethod("parallel", "FraseRSettings", function(object) {
    return(slot(object, "parallel"))
})

#' @export
#' @rdname parallel
setReplaceMethod("parallel", "FraseRSettings", function(object, value) {
    slot(object, "parallel") <- value
    return(object)
})


#'
#'  Get/Set the statistical method to use for P-value calculation
#'  from the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A character string representing the statistical method
#' @examples
#' settings <- createTestFraseRSettings()
#' method(settings)
#' method(settings) <- "betaBinomial"
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname method
setMethod("method", "FraseRSettings", function(object) {
    return(slot(object, "method"))
})

#' @export
#' @rdname method
setReplaceMethod("method", "FraseRSettings", function(object, value) {
    slot(object, "method") <- value
    return(object)
})



## @export
## @rdname method


#'
#'  Get/Set the output folder from the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A path
#' @examples
#' settings <- createTestFraseRSettings()
#' outputFolder(settings)
#' outputFolder(settings) <- tempdir()
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname outputFolder
setMethod("outputFolder", "FraseRSettings", function(object) {
    return(slot(object, "outputFolder"))
})

#' @export
#' @rdname outputFolder
setReplaceMethod("outputFolder", "FraseRSettings", function(object, value) {
    slot(object, "outputFolder") <- value
    return(object)
})


#' 
#' Get/Set the ScanBamParam object from the FraseRSettings object
#' 
#' @param object A FraseRSettings object.
#' @return A ScanBamParam object
#' @examples
#' settings <- createTestFraseRSettings()
#' scanBamParam(settings)
#' scanBamParam(settings) <- ScanBamParam()
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname scanBamParam
setMethod("scanBamParam", "FraseRSettings", function(object) {
    return(slot(object, "bamParams"))
})

#' @export
#' @rdname scanBamParam
setReplaceMethod("scanBamParam", "FraseRSettings", function(object, value) {
    slot(object, "bamParams") <- value
    return(object)
})

#' 
#' Subsetting by indices
#'
#' Providing subsetting by indices through the single-bracket operator
#'
#' @param x A \code{FraseRSettings} object
#' @param i A integer vector
#' @return A subsetted \code{FraseRSettings} object
setMethod("[", c("FraseRSettings", "numeric"), function(x, i) {
    newx <- new("FraseRSettings", 
        sampleData     = sampleData(x)[i],
        bamParams      = scanBamParam(x),
        parallel       = parallel(x),
        method         = method(x),
        strandSpecific = slot(x, "strandSpecific")
    )
    return(newx)
})
setMethod("[", c("FraseRSettings", "logical"), function(x, i) {
    newx <- new("FraseRSettings", 
                sampleData     = sampleData(x)[i],
                bamParams      = scanBamParam(x),
                parallel       = parallel(x),
                method         = method(x),
                strandSpecific = slot(x, "strandSpecific")
    )
    return(newx)
})
