#'
#' @author Christian Mertes
#'
#' This file contains all generic functions
#' some of them are also accessor/getter functions for the FraseRDataSet object

# samples (get the sampleID) functions
#' @export
setGeneric("samples",           function(object) standardGeneric("samples"))
#' @export
setGeneric("samples<-",         signature = "object", function(object, value) standardGeneric("samples<-"))

# condition (get the group/condition of a sample) functions
#' @export
setGeneric("condition",         function(object) standardGeneric("condition"))
#' @export
setGeneric("condition<-",       signature = "object", function(object, value) standardGeneric("condition<-"))

# bamFile (get the bamFile of a sample) functions
#' @export
setGeneric("bamFile",           function(object) standardGeneric("bamFile"))
#' @export
setGeneric("bamFile<-",         signature = "object", function(object, value) standardGeneric("bamFile<-"))

#' @export
setGeneric("name",              function(object) standardGeneric("name"))
#' @export
setGeneric("name<-",            signature = "object", function(object, value) standardGeneric("name<-"))

#' @export
setGeneric("method",            function(object) standardGeneric("method"))
#' @export
setGeneric("method<-",          signature = "object", function(object, value) standardGeneric("method<-"))

#' @export
setGeneric("strandSpecific",    function(object) standardGeneric("strandSpecific"))
#' @export
setGeneric("strandSpecific<-",  signature = "object", function(object, value) standardGeneric("strandSpecific<-"))

# working directory functions
#' @export
setGeneric("workingDir",        function(object) standardGeneric("workingDir"))
#' @export
setGeneric("workingDir<-",      signature = "object", function(object, value) standardGeneric("workingDir<-"))

# scanBamParam functions
#' @export
setGeneric("scanBamParam",      function(object) standardGeneric("scanBamParam"))
#' @export
setGeneric("scanBamParam<-",    signature = "object", function(object, value) standardGeneric("scanBamParam<-"))

# parallel functions
#' @export
setGeneric("parallel",          function(object) standardGeneric("parallel"))
#' @export
setGeneric("parallel<-",        signature = "object", function(object, value) standardGeneric("parallel<-"))

# non spliced reads function
#' @export
setGeneric("nonSplicedReads",   function(object) standardGeneric("nonSplicedReads"))
#' @export
setGeneric("nonSplicedReads<-", signature = "object", function(object, value) standardGeneric("nonSplicedReads<-"))


#'
#'  Getter/Setter for the sampleIDs
#'
#' @param object A FraseRDataSet object.
#' @return A \code{vector} with all sample IDs
#' @examples
#' fds <- createTestFraseRSettings()
#' samples(fds)
#' samples(fds) <- 1:length(fds)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname samples
setMethod("samples", "FraseRDataSet", function(object) {
    return(colData(object)[,"sampleID"])
})

#' @export
#' @rdname samples
setReplaceMethod("samples", "FraseRDataSet", function(object, value) {
    colData(object)[,"sampleID"] <- as.character(value)
    validObject(object)
    return(object)
})


#'
#'  Get the condition for each sample
#'
#' @param object A FraseRDataSet object.
#' @return A \code{vector} with the condition per sample
#'     If the condition column is not set sequence of numbers is returned.
#' @examples
#' fds <- createTestFraseRSettings()
#' condition(fds)
#' condition(fds) <- 1:length(fds)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname condition
setMethod("condition", "FraseRDataSet", function(object) {
    if("condition" %in% colnames(colData(object))){
        return(colData(object)[,"condition"])
    }
    return(1:dim(object)[2])
})

#' @export
#' @rdname condition
setReplaceMethod("condition", "FraseRDataSet", function(object, value) {
    colData(object)[,"condition"] <- value
    validObject(object)
    return(object)
})


#'
#'  Get/Set the bamFile
#'
#' @param object A FraseRDataSet object.
#' @return A \code{vector} with the bamFiles for each sample
#' @examples
#' settings <- createTestFraseRSettings()
#' bamFiles(settings)
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname bamFiles
setMethod("bamFile", "FraseRDataSet", function(object) {
    bamFile <- colData(object)[,"bamFile"]
    if(all(sapply(bamFile, class) == "BamFile")){
        bamFile <- sapply(bamFile, path)
    }
    return(bamFile)
})

#' @export
#' @rdname bamFiles
setReplaceMethod("bamFile", "FraseRDataSet", function(object, value) {
    colData(object)[,"bamFile"] <- value
    validObject(object)
    return(object)
})


#'
#'  Get/Set the parallel param object from the FraseRDataSet object
#'
#' @param object A FraseRDataSet object.
#' @return A parallel param object
#' @examples
#' settings <- createTestFraseRSettings()
#' parallel(settings)
#' parallel(settings) <- SerialParam()
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname parallel
setMethod("parallel", "FraseRDataSet", function(object) {
    return(slot(object, "parallel"))
})

#' @export
#' @rdname parallel
setReplaceMethod("parallel", "FraseRDataSet", function(object, value) {
    slot(object, "parallel") <- value
    validObject(object)
    return(object)
})


#'
#'  Get/Set the name of the analysis
#'
#' @param object A FraseRDataSet object.
#' @return A character string representing the name of the analysis
#' @examples
#' settings <- createTestFraseRSettings()
#' name(settings)
#' name(settings) <- "My Analysis"
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname name
setMethod("name", "FraseRDataSet", function(object) {
    return(slot(object, "name"))
})

#' @export
#' @rdname name
setReplaceMethod("name", "FraseRDataSet", function(object, value) {
    slot(object, "name") <- value
    validObject(object)
    return(object)
})



#'
#'  Get/Set the statistical method to use for P-value calculation
#'  from the FraseRDataSet object
#'
#' @param object A FraseRDataSet object.
#' @return A character string representing the statistical method
#' @examples
#' settings <- createTestFraseRSettings()
#' method(settings)
#' method(settings) <- "betaBinomial"
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname method
setMethod("method", "FraseRDataSet", function(object) {
    return(slot(object, "method"))
})

#' @export
#' @rdname method
setReplaceMethod("method", "FraseRDataSet", function(object, value) {
    slot(object, "method") <- value
    validObject(object)
    return(object)
})


#'
#'  Get/Set the working directory from the FraseRDataSet object
#'
#' @param object A FraseRDataSet object.
#' @return A path
#' @examples
#' settings <- createTestFraseRSettings()
#' workingDir(settings)
#' workingDir(settings) <- tempdir()
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname workingDir
setMethod("workingDir", "FraseRDataSet", function(object) {
    return(slot(object, "workingDir"))
})

#' @export
#' @rdname workingDir
setReplaceMethod("workingDir", "FraseRDataSet", function(object, value) {
    slot(object, "workingDir") <- value
    validObject(object)
    return(object)
})


#'
#'  Get/Set if the analysis is strand specific or not
#'
#' @param object A FraseRDataSet object.
#' @return A logical value if the analysis is strand specific
#' @examples
#' settings <- createTestFraseRSettings()
#' strandSpecific(settings)
#' strandSpecific(settings) <- "betaBinomial"
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname strandSpecific
setMethod("strandSpecific", "FraseRDataSet", function(object) {
    return(slot(object, "strandSpecific"))
})

#' @export
#' @rdname strandSpecific
setReplaceMethod("strandSpecific", "FraseRDataSet", function(object, value) {
    slot(object, "strandSpecific") <- value
    validObject(object)
    return(object)
})


#'
#' Get/Set the ScanBamParam object from the FraseRDataSet object
#'
#' @param object A FraseRDataSet object.
#' @return A ScanBamParam object
#' @examples
#' settings <- createTestFraseRSettings()
#' scanBamParam(settings)
#' scanBamParam(settings) <- ScanBamParam()
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname scanBamParam
setMethod("scanBamParam", "FraseRDataSet", function(object) {
    return(slot(object, "bamParam"))
})

#' @export
#' @rdname scanBamParam
setReplaceMethod("scanBamParam", "FraseRDataSet", function(object, value) {
    slot(object, "bamParam") <- value
    validObject(object)
    return(object)
})


#'
#' accessor for the non spliced reads object within the FraseRDataSet object
#' @export
#' @rdname nonSplicedReads
setMethod("nonSplicedReads", "FraseRDataSet", function(object){
    return(slot(object, "nonSplicedReads"))
})

#'
#' setter for the non spliced reads object within the FraseRDataSet object
#' @export
#' @rdname nonSplicedReads
setReplaceMethod("nonSplicedReads", "FraseRDataSet", function(object, value){
    slot(object, "nonSplicedReads") <- value
    validObject(object)
    return(object)
})

#'
#' Subsetting by indices for junctions
#'
#' Providing subsetting by indices through the single-bracket operator
#'
#' @param x A \code{FraseRDataSet} object
#' @param i A integer vector
#' @return A subsetted \code{FraseRDataSet} object
setMethod("[", c("FraseRDataSet", "ANY", "ANY"), function(x, i, j) {
    if(missing(i) && missing(j)){
        return(x)
    }
    if(missing(i)){
        i <- TRUE
    }
    if(missing(j)){
        j <- TRUE
    }

    # subset the inheritate SE object
    x <- callNextMethod()

    # subset non spliced reads
    nonSplicedReads <- slot(x, "nonSplicedReads")
    if(length(nonSplicedReads) > 0){
        # get selected splice site IDs
        selectedIDs <- unique(unlist(rowData(x)[c("startID", "endID")]))

        # get the selection vector
        iNSR <- rowData(nonSplicedReads)['spliceSiteID'] %in% selectedIDs

        # subset it
        nonSplicedReads <- nonSplicedReads[unlist(iNSR),j]
    }

    newx <- new("FraseRDataSet",
                x,
                name            = name(x),
                method          = method(x),
                parallel        = parallel(x),
                bamParam        = scanBamParam(x),
                strandSpecific  = strandSpecific(x),
                workingDir      = workingDir(x),
                nonSplicedReads = nonSplicedReads
    )
    validObject(newx)
    return(newx)
})


#'
#' Returns the assayNames of FraseR
#'
setMethod("assayNames", "FraseRDataSet", function(x) {
    ans1 <- assayNames(as(x, "SummarizedExperiment"))
    ans2 <- assayNames(nonSplicedReads(x))
    return(c(ans1, ans2))
})


#'
#' Returns the assay corrensonding to the given name/index of the FraseRDataSet
#'
setMethod("assays", "FraseRDataSet", function(x,...){
    return(c(
        callNextMethod(),
        assays(nonSplicedReads(x))
    ))
})
FraseRDataSet.assays.replace <-
    function(x, ..., type="", withDimnames=TRUE, value){
        if(any(names(value) == "")) stop("Name of an assay can not be empty!")
        if(any(duplicated(names(value)))) stop("Assay names need to be unique!")

        # make sure all slots are HDF5
        for(i in seq_along(value)){
            if(!class(value[[i]]) %in% c("HDF5Matrix", "DelayedMatrix")){
                aname <- names(value)[i]
                h5obj <- saveAsHDF5(x, aname, object=value[[i]])
                value[[i]] <- h5obj
            }
        }

        # first replace all existing slots
        n <- length(assayNames(x))
        nj <- n - length(assays(nonSplicedReads(x)))
        ns <- n - nj
        jslots <- value[1:nj]
        sslots <- value[(nj+1):n]
        if(length(value) > n){
            if(!type %in% c("j", "ss")){
                stop(paste("Please set the 'type' option to ",
                        "'j' (Junction) or 'ss' (splice site)."
                ))
            }
            jslots <- c(jslots, value[(1+n):length(value)][type=="j"])
            sslots <- c(sslots, value[(1+n):length(value)][type=="ss"])
        }

        # assign new assays
        value <- jslots
        x <- callNextMethod()
        assays(nonSplicedReads(x), ..., withDimnames=withDimnames) <- sslots

        # validate and return
        validObject(x)
        return(x)
}
setReplaceMethod("assays", c("FraseRDataSet", "SimpleList"),
        FraseRDataSet.assays.replace
)
setReplaceMethod("assays", c("FraseRDataSet", "list"),
        FraseRDataSet.assays.replace
)




