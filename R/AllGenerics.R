
#' @export
setGeneric("samples",           function(object) standardGeneric("samples"))

#' @export
setGeneric("samples<-",         signature = "object", function(object, value) standardGeneric("samples<-"))

#' @export
setGeneric("condition",         function(object) standardGeneric("condition"))

#' @export
setGeneric("condition<-",       signature = "object", function(object, value) standardGeneric("condition<-"))

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

#' @export
setGeneric("workingDir",        function(object) standardGeneric("workingDir"))

#' @export
setGeneric("workingDir<-",      signature = "object", function(object, value) standardGeneric("workingDir<-"))

#' @export
setGeneric("scanBamParam",      function(object) standardGeneric("scanBamParam"))

#' @export
setGeneric("scanBamParam<-",    signature = "object", function(object, value) standardGeneric("scanBamParam<-"))

#' @export
setGeneric("parallel",          function(object) standardGeneric("parallel"))

#' @export
setGeneric("parallel<-",        signature = "object", function(object, value) standardGeneric("parallel<-"))

#' @export
setGeneric("nonSplicedReads",   function(object) standardGeneric("nonSplicedReads"))

#' @export
setGeneric("nonSplicedReads<-", signature = "object", function(object, value) standardGeneric("nonSplicedReads<-"))

#' @rdname results
#' @export
setGeneric("results",   function(x, ...) standardGeneric("results"))


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
#'
#' @rdname samples
#' @export
setMethod("samples", "FraseRDataSet", function(object) {
    return(colData(object)[,"sampleID"])
})

#' @rdname samples
#' @export
setReplaceMethod("samples", "FraseRDataSet", function(object, value) {
    colData(object)[,"sampleID"] <- as.character(value)
    rownames(colData(object)) <- colData(object)[,"sampleID"]
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
#' @param by a character (j or ss) definig if we subset by
#'             junctions or splice sites
#' @return A subsetted \code{FraseRDataSet} object
#' @examples
#'     fds <- countRNAData(createTestFraseRSettings())
#'     fds[1:10,1:10]
#'     fds[,samples(fds) %in% c("sample1", "sample2")]
#'     fds[1:10,by="ss"]
#'
#' @rdname subset
subsetFraseR <- function(x, i, j, by){
    if(missing(i) && missing(j)){
        return(x)
    }
    if(missing(i)){
        i <- TRUE
    }
    if(missing(j)){
        j <- TRUE
    }

    if(is(i, "RangedSummarizedExperiment") | is(i, "GRanges")){
        if(by=="ss"){
            i <- unique(to(findOverlaps(i, nonSplicedReads(x))))
        } else {
            i <- unique(to(findOverlaps(i, x)))
        }
    }

    nsrObj <- nonSplicedReads(x)
    if(length(nsrObj) == 0 & by=="ss"){
        stop("Cannot subset by splice sites, if you not counted them!")
    }
    if(by == "ss"){
        ssIdx <- unlist(rowData(nonSplicedReads(x)[i])["spliceSiteID"])
        i <- as.vector(unlist(rowData(x, typ="psi3")["startID"]) %in% ssIdx |
            unlist(rowData(x, typ="psi3")["endID"]) %in% ssIdx)
    }

    # first subset non spliced reads if needed
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
setMethod("[", c("FraseRDataSet", "ANY", "ANY"),
    function(x, i, j, by=c("j", "ss")) {
        by <- match.arg(by)
        subsetFraseR(x, i, j, by)}
)


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
        assays(nonSplicedReads(x), ..., withDimnames=withDimnames)
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
            side=c("ofInterest", "otherSide"), normalized=FALSE){
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
#' setter for count data
#'
setReplaceMethod("counts", "FraseRDataSet", function(object, type=NULL,
                    side=c("ofInterest", "otherSide"), value){
    side <- match.arg(side)

    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
    } else {
        type <- unlist(
                regmatches(type, gregexpr("psi(3|5|Site)", type, perl=TRUE)))
        aname <- paste0("rawOtherCounts_", type)
    }
    assays(object)[[aname]] <- as.matrix(value)
    validObject(value)
    return(object)
})


options("FraseR-hdf5-chunks"=2000)
options("FraseR-hdf5-cores"=8)
options("FraseR-hdf5-num-chunks"=6)
setAs("DelayedMatrix", "data.table", function(from){
    mc.cores <- min(options()$`FraseR-hdf5-cors`, max(1, detectCores() - 1))
    chunks <- chunk(1:dim(from)[1], options()$`FraseR-hdf5-chunks`)
    fun <- function(x, from) as.data.table(from[x,])

    if(length(chunks) < options()$`FraseR-hdf5-num-chunks`){
        ans <- lapply(chunks, fun, from=from)
    } else {
        ans <- mclapply(chunks, fun, from=from, mc.cores=mc.cores)
        isDT <- sapply(ans, is.data.table)
        if(sum(isDT) != length(chunks)){
            # error happend during extractions, redo
            ans[!isDT] <- lapply(chunks[!isDT], fun, from=from)
        }
    }

    ans <- rbindlist(ans)
    ans
})

setAs("DelayedMatrix", "matrix", function(from){
    as.matrix(as(from, "data.table"))
})

setAs("DataFrame", "data.table", function(from){
    as.data.table(from)
})

setAs("DataFrame", "matrix", function(from){
    as.matrix(as(from, "data.table"))
})

#'
#' retrieve a single sample result object
#' @noRd
resultsSingleSample <- function(sampleID, grs, pvals, zscores, psivals,
                    deltaPsiVals, psiType, fdrCut, zscoreCut, dPsiCut){

    goodCut <- na2false(zscores[,abs(get(sampleID)) >= zscoreCut])
    goodCut <- goodCut & na2false(deltaPsiVals[,abs(get(sampleID)) >= dPsiCut])
    pval    <- pvals[,get(sampleID)]
    p.adj   <- rep(as.numeric(NA), length(pval))
    p.adj[goodCut] <- p.adjust(pval[goodCut], method="hochberg")
    goodCut <- na2false(!is.na(p.adj) & p.adj <= fdrCut & goodCut)

    ans <- granges(grs[goodCut])

    if(!any(goodCut)){
        return(ans)
    }

    # extract data
    mcols(ans)$type        <- psiType
    mcols(ans)$sampleID    <- sampleID
    mcols(ans)$hgnc_symbol <- mcols(grs[goodCut])$hgnc_symbol
    mcols(ans)$zscore      <- round(zscores[goodCut,get(sampleID)], 3)
    mcols(ans)$psiValue    <- round(psivals[goodCut,get(sampleID)], 3)
    mcols(ans)$pvalue      <- pval[goodCut]
    mcols(ans)$deltaPsi    <- round(deltaPsiVals[goodCut,get(sampleID)], 3)

    # correct for multiple testing
    p <- mcols(ans)$pvalue
    n <- length(grs)
    mcols(ans)$p.adj       <- p.adj[goodCut]

    return(ans[order(mcols(ans)$pvalue)])
}

FraseR.results <- function(x, sampleIDs, fdrCut, zscoreCut, dPsiCut,
                    psiType, redo){

    # check input
    stopifnot(is.numeric(fdrCut)    && fdrCut    <= 1   && fdrCut    >= 0)
    stopifnot(is.numeric(dPsiCut)   && dPsiCut   <= 1   && dPsiCut   >= 0)
    stopifnot(is.numeric(zscoreCut) && zscoreCut <= 100 && zscoreCut >= 0)

    stopifnot(is(x, "FraseRDataSet"))
    stopifnot(all(sampleIDs %in% samples(x)))

    # check if we extacted it already
    fdrCut    <- round(fdrCut, 3)
    dPsiCut   <- round(dPsiCut, 2)
    zscoreCut <- round(zscoreCut, 2)

    # get the name of the result
    if(length(sampleIDs)==1){
        sampleStr <- sampleIDs
    } else if(all(samples(x) %in% sampleIDs)){
        sampleStr <- "all"
    } else {
        sampleStr <- paste(sort(sampleIDs), collapse=", ")
    }
    resName <- sprintf(paste0("Results for samples: %s, type: %s, ",
                "and cutoffs: padjc <= %g & abs(zc) >= %g & abs(dpsic) >= %g"),
        sampleStr, paste(psiType, collapse=" + "), fdrCut, zscoreCut, dPsiCut
    )

    # return cache if there
    if(redo==FALSE & resName %in% names(metadata(x))){
        message(date(), ": Used result cache: ", resName)
        return(metadata(x)[[resName]])
    }

    resultsls <- sapply(psiType, function(type){
        message(date(), ": Collecting results for: ", type)

        tested <- mcols(x, type=type)[[paste0(type, "_tested")]]
        tested <- na2false(tested)

        # get ranges
        grs <- rowRanges(if(type=="psiSite") nonSplicedReads(x) else x)
        grs <- grs[tested]


        # extract values
        pvals <- as(Class="data.table",
                object=assays(x)[[paste0("pvalue_", type)]][tested,])
        zscores <- as(Class="data.table",
                object=assays(x)[[paste0("zscore_", type)]][tested,])
        psivals <- as(Class="data.table", object=assays(x)[[type]][tested,])
        deltaPsiVals <- as(Class="data.table", object=assays(x)[[
                paste0("delta_", type)]][tested,])

        # create reult table
        sampleRes <- sapply(sampleIDs,
                resultsSingleSample, grs=grs, pvals=pvals,
                zscores=zscores, psiType=type, psivals=psivals,
                deltaPsiVals=deltaPsiVals,
                fdrCut=fdrCut, zscoreCut=zscoreCut, dPsiCut=dPsiCut)

        # return combined result
        return(unlist(GRangesList(sampleRes)))
    })

    # merge results
    ans <- unlist(GRangesList(resultsls))

    # save it for later into the orignial object
    metadata(x)[[resName]] <- ans

    # try to set it in the user environemnt
    try({
        sysCall <- sys.calls()
        resCall <- deparse(sysCall[[length(sysCall)-2]])
        fdsName <- strsplit(resCall, "[(,]")[[1]][2]
        if(fdsName %in% ls(envir = parent.frame(3))){
            message(date(), ": Caching the results in the object: '", fdsName,
                    "' under the metadata with the name: '", resName, "'.")
            assign(fdsName, x, envir = parent.frame(3))
        }
    })

    # return only the results
    return(ans)
}

#'
#' Extracting results
#'
#' The result function extracts the results from the given analysis object
#' based on the given options and cutoffs.
#'
#' @rdname results
#' @export
setMethod("results", "FraseRDataSet", function(x, sampleIDs=samples(x),
                    fdrCut=0.05, zscoreCut=2, redo=FALSE,
                    dPsiCut=0.01, psiType=c("psi3", "psi5", "psiSite")){
    FraseR.results(x, sampleIDs=sampleIDs, fdrCut=fdrCut, zscoreCut=zscoreCut,
            redo=redo, dPsiCut=dPsiCut,
            psiType=match.arg(psiType, several.ok=TRUE))
})

#'
#' Mapping of chromosome names
#'
#' @examples
#'
#' fds <- countRNAData(createTestFraseRSettings())
#'
#' seqlevels(mapSeqlevels(fds, style="UCSC"))
#' seqlevels(mapSeqlevels(fds, style="Ensembl"))
#' seqlevels(mapSeqlevels(fds, style="dbSNP"))
#'
#' @export
mapSeqlevels <- function(fds, style="UCSC", ...){

    mappings <- na.omit(GenomeInfoDb::mapSeqlevels(seqlevels(fds), style, ...))

    if(length(mappings) != length(seqlevels(fds))){
        message(date(), ": Drop non standard chromosomes for compatibility.")
        fds <- keepStandardChromosomes(fds)
        nonSplicedReads(fds) <- keepStandardChromosomes(nonSplicedReads(fds))
        validObject(fds)
    }
    fds <- fds[as.vector(seqnames(fds)) %in% names(mappings)]

    seqlevels(fds) <- as.vector(mappings)
    seqlevels(nonSplicedReads(fds)) <- as.vector(mappings)

    return(fds)
}


