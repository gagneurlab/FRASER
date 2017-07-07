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
#' samples(fds) <- 1:dim(fds)[2]
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
#' condition(fds) <- 1:dim(fds)[2]
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
#' @return A \code{vector} with the bamFile for each sample
#' @examples
#' settings <- createTestFraseRSettings()
#' bamFile(settings)
#' bamFile(settings) <- file.path("bamfiles", samples(settings), "rna-seq.bam")
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @rdname bamFile
setMethod("bamFile", "FraseRDataSet", function(object) {
    bamFile <- colData(object)[,"bamFile"]
    if(all(sapply(bamFile, class) == "BamFile")){
        bamFile <- sapply(bamFile, path)
    }
    return(bamFile)
})

#' @export
#' @rdname bamFile
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
#' method(settings) <- "betaBin"
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
#' strandSpecific(settings) <- TRUE
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
#' scanBamParam(settings) <- ScanBamParam(mapqFilter=30)
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
#' @param i A integer vector to subset the rows/ranges
#' @param j A integer vector to subset the columns/samples
#' @return A subsetted \code{FraseRDataSet} object
#' @examples
#'     fds <- countRNAData(createTestFraseRSettings())
#'     fds[1:10,1:10]
#'     fds[,samples(fds) %in% c("sample1", "sample2")]
#'
#' @rdname subset
subsetFraseR <- function(x, i, j){
    if(missing(i) && missing(j)){
        return(x)
    }
    if(missing(i)){
        i <- TRUE
    }
    if(missing(j)){
        j <- TRUE
    }

    # first subset non spliced reads if needed
    nsrObj <- nonSplicedReads(x)
    if(length(nsrObj) > 0){
        # get selected splice site IDs
        selectedIDs <- unique(unlist(
            rowData(x, type="j")[i, c("startID", "endID")]
        ))

        # get the selection vector
        idxNSR <- rowData(x, type="ss")[['spliceSiteID']] %in% selectedIDs

        # subset it
        nsrObj <- nsrObj[idxNSR,j]
    }

    # subset the inheritate SE object
    if(length(x) == 0){
        i <- NULL
    }
    subX <- as(as(x, "RangedSummarizedExperiment")[i,j], "FraseRDataSet")

    # create new FraseRDataSet object
    newx <- new("FraseRDataSet",
                subX,
                name            = name(x),
                method          = method(x),
                parallel        = parallel(x),
                bamParam        = scanBamParam(x),
                strandSpecific  = strandSpecific(x),
                workingDir      = workingDir(x),
                nonSplicedReads = nsrObj
    )
    validObject(newx)
    return(newx)
}
setMethod("[", c("FraseRDataSet", "ANY", "ANY"), function(x, i, j) {
    subsetFraseR(x,i,j)
})


#'
#' Returns the assayNames of FraseR
#'
setMethod("assayNames", "FraseRDataSet", function(x) {
    return(c(
        assayNames(as(x, "SummarizedExperiment")),
        assayNames(nonSplicedReads(x))
    ))
})


#'
#' Returns the assay corrensonding to the given name/index of the FraseRDataSet
#'
setMethod("assays", "FraseRDataSet", function(x, ..., type=NULL, withDimnames=TRUE){
    return(c(
        callNextMethod(),
        assays(nonSplicedReads(x, ...))
    ))
})
FraseRDataSet.assays.replace <-
            function(x, ..., type=NULL, withDimnames=TRUE, value){
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
    nj <- names(value) %in% assayNames(as(x, "SummarizedExperiment"))
    ns <- names(value) %in% assayNames(nonSplicedReads(x))
    jslots <- value[nj]
    sslots <- value[ns]

    # add new slots if there are some
    if(sum(!(nj | ns)) > 0){
        type <- sapply(type, checkReadType, fds=x)
        jslots <- c(jslots, value[!(nj | ns)][type=="j"])
        sslots <- c(sslots, value[!(nj | ns)][type=="ss"])
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
setReplaceMethod("assays", c("FraseRDataSet", "DelayedMatrix"),
        FraseRDataSet.assays.replace
)

#'
#' retrive the length of the object (aka number of junctions)
#'
setMethod("length", "FraseRDataSet", function(x) callNextMethod())

#'
#' getter and setter for mcols
#'
setMethod("mcols", "FraseRDataSet", function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j")  return(callNextMethod())
    if(type=="ss") return(mcols(nonSplicedReads(x), ...))
})
setReplaceMethod("mcols", "FraseRDataSet", function(x, type=NULL, ..., value){
    type <- checkReadType(x, type)
    if(type=="j") return(callNextMethod())
    if(type=="ss"){
        mcols(nonSplicedReads(x), ...) <- value
        return(x)
    }
})


#'
#' getter for count data
#'
setMethod("counts", "FraseRDataSet", function(object, type=NULL,
            side=c("ofInterest", "otherSide")){
    side <- match.arg(side)
    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
        if(!aname %in% assayNames(object)){
            stop("Missing rawCounts. Please count your data first. ",
                 "And then try again."
            )
        }
        return(assays(object)[[aname]])
    }

    # extract psi value from type
    type <- unlist(regmatches(type, gregexpr("psi(3|5|Site)", type, perl=TRUE)))
    aname <- paste0("rawOtherCounts_", type)
    if(!aname %in% assayNames(object)){
        stop("Missing rawOtherCounts. Please calculate PSIValues first. ",
             "And then try again."
        )
    }
    return(assays(object)[[aname]])
})

#'
#' convertion of Delayed Matrix objects into a data.table
#'
setAs("DelayedMatrix", "data.table", function(from){
    mc.cores=min(24, max(1, detectCores() - 1))
    perChunk=5000
    chunks <- chunk(1:dim(from)[1], perChunk)
    ans <- mclapply(chunks, mc.cores=mc.cores,
        FUN=function(x) as.data.table(from[x,])
    )
    ans <- rbindlist(ans)
    ans
})

#'
#' convertion of Delayed Matrix objects into a matrix
#'
setAs("DelayedMatrix", "matrix", function(from){
    #mc.cores=min(24, max(1, detectCores() - 1))
    #perChunk=30000
    #chunks <- chunk(1:dim(from)[1], perChunk)
    #ans <- mclapply(chunks, mc.cores=mc.cores,
    #                FUN=function(x) as.matrix(from[x,])
    #)
    #ans1 <- rbindlist(ans1)
    #ans
    as.matrix(as(from, "data.table"))
})

#'
#' retrieve a single sample result object
#' @noRd
resultsSingleSample <- function(sampleID, grs, pvals, zscores, psivals,
                                psiType, pvalueCut, zscoreCut){
    goodCut <- na2false(
        pvals[,get(sampleID) <= pvalueCut] &
            zscores[,abs(get(sampleID)) >= zscoreCut]
    )

    ans <- granges(grs[goodCut])

    if(!any(goodCut)){
        return(ans)
    }

    # extract data
    mcols(ans)$type        <- psiType
    mcols(ans)$sampleID    <- sampleID
    mcols(ans)$hgnc_symbol <- mcols(grs[goodCut])$hgnc_symbol
    mcols(ans)$zscore      <- round(zscores[goodCut,get(sampleID)], 2)
    mcols(ans)$psiValue    <- round(psivals[goodCut,get(sampleID)], 2)
    mcols(ans)$pvalue      <- pvals[goodCut,get(sampleID)]

    # correct for multiple testing
    p <- mcols(ans)$pvalue
    n <- length(goodCut)
    mcols(ans)$p.hochberg  <- p.adjust(p, n=n, method="hochberg")
    mcols(ans)$fdr         <- p.adjust(p, n=n, method="fdr")

    return(ans[order(mcols(ans)$pvalue)])
}

#'
#' obtain the results for the given analysis pipeline
#' @export
results <- function(fds, sampleIDs=samples(fds), pvalueCut=1e-5, zscoreCut=2,
                    psiType=c("psi3", "psi5", "psiSite"), redo=FALSE){

    # check input
    stopifnot(is.numeric(pvalueCut) && pvalueCut <= 1 && pvalueCut >= 0)
    stopifnot(is.numeric(zscoreCut) && zscoreCut <= 100 && zscoreCut >= 0)
    stopifnot(is(fds, "FraseRDataSet"))
    stopifnot(all(sampleIDs %in% samples(fds)))
    psiType <- match.arg(psiType, several.ok=TRUE)

    # check if we extacted it already
    pvalueCut <- round(pvalueCut, 8)
    zscoreCut <- round(zscoreCut, 2)

    # get the name of the result
    if(length(sampleIDs)==1){
        sampleStr <- sampleIDs
    } else if(length(sampleIDs == dim(fds)[2])) {
        sampleStr <- "all"
    } else {
        sampleStr <- deparse(which(samples(fds) %in% sampleIDs))
    }
    resName <- sprintf(
        "Results for samples: %s and cutoffs: pc <= %g & abs(zc) >= %g",
        sampleStr, pvalueCut, zscoreCut
    )

    # return cache if there
    if(redo==FALSE & resName %in% names(metadata(fds))){
        message(date(), ": Used result cache: ", resName)
        return(metadata(fds)[[resName]])
    }

    resultsls <- sapply(psiType, function(type){
        message(date(), ": Collecting results for: ", type)

        tested <- mcols(fds, type=type)[[paste0(type, "_tested")]]
        tested <- na2false(tested)
        pvals <- as(Class="data.table",
                object=assays(fds)[[paste0("pvalue_", type)]][tested,]
        )

        zscores <- as(Class="data.table",
                object=assays(fds)[[paste0("zscore_", type)]][tested,]
        )

        psivals <- as(Class="data.table", object=assays(fds)[[type]][tested,])

        grs <- rowRanges(if(type=="psiSite") nonSplicedReads(fds) else fds)
        grs <- grs[tested]

        sampleRes <- sapply(sampleIDs,
               resultsSingleSample, grs=grs, pvals=pvals,
               zscores=zscores, psiType=type, psivals=psivals,
               pvalueCut=pvalueCut, zscoreCut=zscoreCut
        )

        results <- unlist(GRangesList(sampleRes))
    })

    ans <- unlist(GRangesList(resultsls))

    # save it for later into the orignial object
    fdsObjName <- deparse(substitute(fds))
    eval(parse(text = sprintf("%s <<- %s",
            paste0("metadata(", fdsObjName, ")[['",resName, "']]"),
            "ans"
    )))

    # return only the results
    return(ans)
}


