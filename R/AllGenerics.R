asSE <- function(x){
    if(is(x, "RangedSummarizedExperiment")){
        return(as(x, "RangedSummarizedExperiment"))
    }
    as(x, "SummarizedExperiment")
}

asFDS <- function(x){
    return(as(x, "FraseRDataSet"))
}


#' @rdname samples
#' @export
setGeneric("samples",           
            function(object) standardGeneric("samples"))

#' @rdname samples
#' @export
setGeneric("samples<-",         signature = "object", 
            function(object, value) standardGeneric("samples<-"))

#' @rdname condition
#' @export
setGeneric("condition",         
            function(object) standardGeneric("condition"))

#' @rdname condition
#' @export
setGeneric("condition<-",       signature = "object", 
            function(object, value) standardGeneric("condition<-"))

#' @rdname bamFile
#' @export
setGeneric("bamFile",           
            function(object) standardGeneric("bamFile"))

#' @rdname bamFile
#' @export
setGeneric("bamFile<-",         signature = "object", 
            function(object, value) standardGeneric("bamFile<-"))

#' @rdname name
#' @export
setGeneric("name",              
            function(object) standardGeneric("name"))

#' @rdname name
#' @export
setGeneric("name<-",            signature = "object", 
            function(object, value) standardGeneric("name<-"))

#' @rdname strandSpecific
#' @export
setGeneric("strandSpecific",    
            function(object) standardGeneric("strandSpecific"))

#' @rdname strandSpecific
#' @export
setGeneric("strandSpecific<-",  signature = "object", 
            function(object, value) standardGeneric("strandSpecific<-"))

#' @rdname workingDir
#' @export
setGeneric("workingDir",        
            function(object) standardGeneric("workingDir"))

#' @rdname workingDir
#' @export
setGeneric("workingDir<-",      signature = "object", 
            function(object, value) standardGeneric("workingDir<-"))

#' @rdname scanBamParam
#' @export
setGeneric("scanBamParam",      
            function(object) standardGeneric("scanBamParam"))

#' @rdname scanBamParam
#' @export
setGeneric("scanBamParam<-",    signature = "object", 
            function(object, value) standardGeneric("scanBamParam<-"))

#' @rdname nonSplicedReads
#' @export
setGeneric("nonSplicedReads",   
            function(object) standardGeneric("nonSplicedReads"))

#' @rdname nonSplicedReads
#' @export
setGeneric("nonSplicedReads<-", signature = "object", 
            function(object, value) standardGeneric("nonSplicedReads<-"))

#' @rdname results
#' @return GRanges object
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))


#'
#'  Getter/Setter for the sampleIDs
#'
#' @param object A FraseRDataSet object.
#' @param value The new sample IDs.
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
#' @param value A vector with the conditions that should be assigned.
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
    return(seq_len(dim(object)[2]))
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
#' @param value The new values (file.paths or objects of class BamFile) for the 
#' bamfiles. 
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
    if(all(vapply(bamFile, is, class2="BamFile", FUN.VALUE=logical(1)))){
        bamFile <- vapply(bamFile, path, "")
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
#'  Get/Set the name of the analysis
#'
#' @param object A FraseRDataSet object.
#' @param value The new name of the analysis.
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
#'  Get/Set the working directory from the FraseRDataSet object
#'
#' @param object A FraseRDataSet object.
#' @param value The new working directory.
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
#' @param value The new value for the strand specificity of the analysis.
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
#' @param value The new ScanBamParam.
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


#' @export
#' @rdname nonSplicedReads
setMethod("nonSplicedReads", "FraseRDataSet", function(object){
    return(slot(object, "nonSplicedReads"))
})

#'
#' Getter/setter for the non spliced reads object within the FraseRDataSet 
#' object
#' 
#' @param object A FraseRDataSet object.
#' @param value A RangedSummarizedExperiment object containing the counts for 
#' the non spliced reads overlapping splice sites in the fds.
#' @return RangedSummarizedExperiment (getter) or FraseRDataSet (setter)
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
subset.FraseR <- function(x, i, j, by=c("j", "ss")){
    if(length(by) == 1){
        by <- whichReadType(x, by)
    }
    by <- match.arg(by)

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
    ssIdx <- NULL
    if(by == "ss"){
        ssIdx <- unlist(rowData(nonSplicedReads(x)[i])["spliceSiteID"])
        i <- as.vector(unlist(rowData(x, typ="psi3")["startID"]) %in% ssIdx |
            unlist(rowData(x, typ="psi3")["endID"]) %in% ssIdx)
    }

    # first subset non spliced reads if needed
    if(length(nsrObj) > 0){
        if(is.null(ssIdx)){
            # get selected splice site IDs
            ssIdx <- unique(unlist(
                rowData(x, type="j")[i, c("startID", "endID")]
            ))
        }

        # get the selection vector
        idxNSR <- rowData(x, type="ss")[['spliceSiteID']] %in% ssIdx

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
            bamParam        = scanBamParam(x),
            strandSpecific  = strandSpecific(x),
            workingDir      = workingDir(x),
            nonSplicedReads = nsrObj
    )
    validObject(newx)
    return(newx)
}
#' @rdname subset
#' @export
setMethod("[", c("FraseRDataSet", "ANY", "ANY"), subset.FraseR)


#'
#' Returns the assayNames of FraseR
#' 
#' @param x FraseRDataSet
#' 
#' @return Character vector
#' @export
setMethod("assayNames", "FraseRDataSet", function(x) {
    return(c(
        assayNames(asSE(x)),
        assayNames(nonSplicedReads(x))
    ))
})


#'
#' Returns the assay corrensonding to the given name/index of the FraseRDataSet
#'
#' @param x FraseRDataSet
#' @param ... Parameters passed on to SummarizedExperiment::assays() 
#' @param withDimnames Passed on to SummarizedExperiment::assays() 
#'
#' @return (Delayed) matrix.
#' @export
setMethod("assays", "FraseRDataSet", function(x, ..., withDimnames=TRUE){
    return(c(
        assays(asSE(x), ..., withDimnames=withDimnames),
        assays(nonSplicedReads(x), ..., withDimnames=withDimnames)
    ))
})
FraseRDataSet.assays.replace <-
            function(x, ..., HDF5=TRUE, type=NULL, withDimnames=TRUE, value){
    if(any(names(value) == "")) stop("Name of an assay can not be empty!")
    if(any(duplicated(names(value)))) stop("Assay names need to be unique!")
    if(is.null(type)){
        type <- names(value)
    }
    # make sure all slots are HDF5
    if(isTRUE(HDF5)){
        for(i in seq_along(value)){
            if(!class(value[[i]]) %in% c("HDF5Matrix", "DelayedMatrix")){
                aname <- names(value)[i]
                h5obj <- saveAsHDF5(x, aname, object=value[[i]])
                value[[i]] <- h5obj
            }
        }
    }

    # first replace all existing slots
    nj <- names(value) %in% assayNames(asSE(x))
    ns <- names(value) %in% assayNames(nonSplicedReads(x))
    jslots <- value[nj]
    sslots <- value[ns]

    # add new slots if there are some
    if(sum(!(nj | ns)) > 0){
        type <- vapply(type, checkReadType, fds=x, FUN.VALUE="")
        jslots <- c(jslots, value[!(nj | ns) & type=="j"])
        sslots <- c(sslots, value[!(nj | ns) & type=="ss"])
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
#' @param x FraseRDataSet
#'
#' @return Length of the object.
#' @export
setMethod("length", "FraseRDataSet", function(x) callNextMethod())

#'
#' getter and setter for mcols
#' 
#' @param x FraseRDataSet
#' @param type The psi type.
#' @param ... Further parameters passed to mcols
#' 
#' @return mcols
#' @examples 
#'   fds <- createTestFraseRDataSet()
#'   
#'   mcols(fds, type="psi5")
#'   mcols(fds, type="psiSite")
#' @export
FraseR.mcols.get <- function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j"){
        return(mcols(asSE(x), ...))
    }
    mcols(nonSplicedReads(x), ...)
}
FraseR.mcols.replace <- function(x, type=NULL, ..., value){
    type <- checkReadType(x, type)
    if(type=="j") {
        return(callNextMethod())
        # se <- asSE(x)
        # mcols(se, ...) <- value
        # return(asFDS(x))
    }
    mcols(nonSplicedReads(x), ...) <- value
    return(x)
}
setMethod("mcols", "FraseRDataSet", FraseR.mcols.get)
setReplaceMethod("mcols", "FraseRDataSet", FraseR.mcols.replace)

#'
#' Getter and setter for rowRanges
#' 
#' @param x FraseRDataSet
#' @param type The psi type.
#' @param ... Parameters passed on to rowRanges()
#' 
#' @return GRanges object
#' @examples 
#'   fds <- createTestFraseRDataSet()
#'   
#'   rowRanges(fds)
#'   rowRanges(fds, type="psiSite")
#' 
#' @return FraseRDataSet
FraseR.rowRanges.get <- function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j")  return(callNextMethod())
    if(type=="ss") return(rowRanges(nonSplicedReads(x), ...))
}
FraseR.rowRanges.replace <- function(x, type=NULL, ..., value){
    type <- checkReadType(x, type)
    if(type=="j") return(callNextMethod())
    rowRanges(nonSplicedReads(x), ...) <- value
    return(x)
}

setMethod("rowRanges", "FraseRDataSet", FraseR.rowRanges.get)
setReplaceMethod("rowRanges", "FraseRDataSet", FraseR.rowRanges.replace)

#'
#' Getter/setter for count data
#'
#' @param fds,object FraseRDataSet
#' @param type The psi type.
#' @param side "ofInterest" for junction counts, "other" for sum of counts of 
#' all other junctions at the same donor site (psi5) or acceptor site (psi3), 
#' respectively. 
#' @param value An integer matrix containing the counts.
#' @param ... Further parameters that are passed to assays(object,...)
#' 
#' @return FraseRDataSet
#' @rdname counts
#' @examples 
#'  fds <- createTestFraseRDataSet()
#'  
#'  counts(fds, type="psi5", side="ofInterest")
#'  counts(fds, type="psi5", side="other")
#'  head(K(fds, type="psi3"))
#'  head(N(fds, type="psi3"))
#'  
setMethod("counts", "FraseRDataSet", function(object, type=NULL,
            side=c("ofInterest", "otherSide")){
    side <- match.arg(side)
    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
        if(!aname %in% assayNames(object)){
            stop(paste0("Missing rawCounts for '", type, "'. ",
                    "Please count your data first. And then try again."))
        }
        return(assays(object)[[aname]])
    }

    # extract psi value from type
    type <- unlist(regmatches(type, gregexpr("psi(3|5|Site)", type, perl=TRUE)))
    if(length(type) == 0){
        stop(paste0("Please provide a correct psi type: psi5, psi3, and ",
                    "psiSite. Not the given one: '", type, "'."))
    }
    aname <- paste0("rawOtherCounts_", type)
    if(!aname %in% assayNames(object)){
        stop(paste0("Missing rawOtherCounts for type '", type, "'.",
                "Please calculate PSIValues first. And then try again."))
    }
    return(assays(object)[[aname]])
})

#'
#' setter for count data
#' 
#' @rdname counts
setReplaceMethod("counts", "FraseRDataSet", function(object, type=NULL,
                    side=c("ofInterest", "otherSide"), ..., value){
    side <- match.arg(side)

    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
    } else {
        type <- unlist(
                regmatches(type, gregexpr("psi(3|5|Site)", type, perl=TRUE)))
        aname <- paste0("rawOtherCounts_", type)
    }
    assays(object, ...)[[aname]] <- value
    validObject(value)
    return(object)
})


setAs("DelayedMatrix", "data.table", function(from){ 
    as.data.table(from) })
setAs("DataFrame",     "data.table", function(from){ 
    as.data.table(from) })
setAs("DelayedMatrix", "matrix", function(from){ 
    as.matrix(as(from, "data.table")) })
setAs("DataFrame", "matrix", function(from){
    as.matrix(as(from, "data.table")) })

#'
#' retrieve a single sample result object
#' @noRd
resultsSingleSample <- function(sampleID, gr, pvals, padjs, zscores, psivals,
                rawCts, rawTotalCts, deltaPsiVals, muPsi, psiType, fdrCut,
                zscoreCut, dPsiCut, rowMeansK, rowMeansN, minCount){

    zscore  <- zscores[,sampleID]
    dpsi    <- deltaPsiVals[,sampleID]
    pval    <- pvals[,sampleID]
    padj    <- padjs[,sampleID]

    goodCut <- !logical(length(zscore))
    if(!is.na(zscoreCut)){
        goodCut <- goodCut & na2default(abs(zscore) >= zscoreCut, TRUE)
    }
    if(!is.na(dPsiCut)){
        goodCut <- goodCut & na2default(abs(dpsi) >= dPsiCut, TRUE)
    }
    if(!is.na(fdrCut)){
        goodCut <- goodCut & na2false(padj <= fdrCut)
    }
    if(!is.na(minCount)){
        goodCut <- goodCut & rawTotalCts[,sampleID] >= minCount
    }

    ans <- granges(gr[goodCut])

    if(!any(goodCut)){
        return(ans)
    }

    if(!"hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(gr)$hgnc_symbol <- NA_character_
    }

    # extract data
    mcols(ans)$sampleID        <- Rle(sampleID)
    mcols(ans)$hgncSymbol      <- Rle(mcols(gr[goodCut])$hgnc_symbol)
    mcols(ans)$type            <- Rle(psiType)
    mcols(ans)$pValue          <- signif(pval[goodCut], 5)
    mcols(ans)$padjust         <- signif(padj[goodCut], 5)
    mcols(ans)$zScore          <- Rle(round(zscore[goodCut], 2))
    mcols(ans)$psiValue        <- Rle(round(psivals[goodCut,sampleID], 2))
    mcols(ans)$deltaPsi        <- Rle(round(dpsi[goodCut], 2))
    mcols(ans)$meanCounts      <- Rle(round(rowMeansK[goodCut], 2))
    mcols(ans)$meanTotalCounts <- Rle(round(rowMeansN[goodCut], 2))
    mcols(ans)$counts          <- Rle(rawCts[goodCut, sampleID])
    mcols(ans)$totalCounts     <- Rle(rawTotalCts[goodCut, sampleID])

    return(ans[order(mcols(ans)$pValue)])
}

FraseR.results <- function(x, sampleIDs, fdrCutoff, zscoreCutoff, dPsiCutoff,
                    psiType, BPPARAM=bpparam(), maxCols=20, minCount){

    # check input
    checkNaAndRange(fdrCutoff,    min=0, max=1,   scalar=TRUE, na.ok=TRUE)
    checkNaAndRange(dPsiCutoff,   min=0, max=1,   scalar=TRUE, na.ok=TRUE)
    checkNaAndRange(zscoreCutoff, min=0, max=100, scalar=TRUE, na.ok=TRUE)
    checkNaAndRange(minCount,     min=0, max=Inf, scalar=TRUE, na.ok=TRUE)

    stopifnot(is(x, "FraseRDataSet"))
    stopifnot(all(sampleIDs %in% samples(x)))

    resultsls <- bplapply(psiType, BPPARAM=BPPARAM, function(type){
        message(date(), ": Collecting results for: ", type)
        currentType(x) <- type
        gr <- rowRanges(x, type=type)

        # first get row means
        rowMeansK <- rowMeans(K(x, type=type))
        rowMeansN <- rowMeans(N(x, type=type))

        # then iterate by chunk
        chunkCols <- getMaxChunks2Read(fds=x, assayName=type, max=maxCols)
        sampleChunks <- getSamplesByChunk(fds=x, sampleIDs=sampleIDs,
                chunkSize=chunkCols)

        ans <- lapply(seq_along(sampleChunks), function(idx){
            message(date(), ": Process chunk: ", idx, " for: ", type)
            sc <- sampleChunks[[idx]]
            tmp_x <- x[,sc]

            # extract values
            rawCts       <- as.matrix(K(tmp_x))
            rawTotalCts  <- as.matrix(N(tmp_x))
            pvals        <- as.matrix(pVals(tmp_x))
            padjs        <- as.matrix(padjVals(tmp_x))
            zscores      <- as.matrix(zScores(tmp_x))
            psivals      <- as.matrix(assay(tmp_x, type))
            muPsi        <- as.matrix(predictedMeans(tmp_x))
            deltaPsiVals <- psivals - muPsi

            if(length(sc) == 1){
                colnames(pvals) <- sc
                colnames(padjs) <- sc
                colnames(zscores) <- sc
                colnames(deltaPsiVals) <- sc
            }
            
            # create reult table
            sampleRes <- lapply(sc,
                    resultsSingleSample, gr=gr, pvals=pvals, padjs=padjs,
                    zscores=zscores, psiType=type, psivals=psivals,
                    deltaPsiVals=deltaPsiVals, muPsi=muPsi, rawCts=rawCts,
                    rawTotalCts=rawTotalCts, fdrCut=fdrCutoff,
                    zscoreCut=zscoreCutoff, dPsiCut=dPsiCutoff,
                    rowMeansK=rowMeansK, rowMeansN=rowMeansN, 
                    minCount=minCount)

            # return combined result
            return(unlist(GRangesList(sampleRes)))
        })

        unlist(GRangesList(ans))
    })

    # merge results
    ans <- unlist(GRangesList(resultsls))

    # sort it if existing
    if(length(ans) > 0){
        ans <- ans[order(ans$pValue)]
    }

    # return only the results
    return(ans)
}

#'
#' Extracting results
#'
#' The result function extracts the results from the given analysis object
#' based on the given options and cutoffs.
#'
#' @param x FraseRDataSet
#' @param sampleIDs A vector of sample IDs for which results should be 
#' retrieved
#' @param padjCutoff The FDR cutoff to be applied.
#' @param zScoreCutoff The z-score cutoff to be applied.
#' @param deltaPsiCutoff The cutoff on delta psi.
#' @param minCount The minimum count value of the total coverage of an intron 
#' to be considered as significant.
#' result
#' @param psiType The psi types for which the results should be retrieved.
#' @param BPPARAM The BiocParallel parameter.
#' @param ... Further parameters that are passed on.
#'
#' @return GRanges object containing significant results.
#' 
#' @rdname results
#' @examples
#' # get data, fit and compute p-values and z-scores
#' fds <- createTestFraseRDataSet()
#' 
#' # extract results: for this example dataset, padjust cutoff of 0.15 used to
#' # get at least one result and show the output
#' res <- results(fds, padjCutoff=0.15, zScoreCutoff=NA, deltaPsiCutoff=0.1)
#' res
#'
#' @export
setMethod("results", "FraseRDataSet", function(x, sampleIDs=samples(x),
                    padjCutoff=0.05, zScoreCutoff=NA, deltaPsiCutoff=0.3,
                    minCount=5, psiType=c("psi3", "psi5", "psiSite"),
                    BPPARAM=bpparam()){
    FraseR.results(x, sampleIDs=sampleIDs, fdrCutoff=padjCutoff,
            zscoreCutoff=zScoreCutoff, dPsiCutoff=deltaPsiCutoff,
            minCount=minCount, psiType=match.arg(psiType, several.ok=TRUE),
            BPPARAM=BPPARAM)
})

resultsByGenes <- function(res, geneColumn="hgncSymbol", method="BY"){
    # sort by pvalue
    res <- res[order(res$pValue)]

    # extract subset
    if(is(res, "GRanges")){
        ans <- as.data.table(mcols(res)[,c(geneColumn, "pValue", "sampleID")])
        colnames(ans) <- c("features", "pval", "sampleID")
    } else {
        ans <- featureNames <- res[,.(
                features=get(geneColumn), pval=pvalue, sampleID=sampleID)]
    }

    # remove NAs
    naIdx <- ans[,is.na(features)]
    ansNoNA <- ans[!is.na(features)]

    # compute pvalues by gene
    ansNoNA[,pByFeature:=min(p.adjust(pval, method="holm")),
            by="sampleID,features"]

    # subset to lowest pvalue by gene
    dupIdx <- duplicated(ansNoNA[,.(features,sampleID)])
    ansGenes <- ansNoNA[!dupIdx]

    # compute FDR
    ansGenes[,fdrByFeature:=p.adjust(pByFeature, method=method),
            by="sampleID"]

    # get final result table
    finalAns <- res[!naIdx][!dupIdx]
    finalAns$pValueGene  <- ansGenes$pByFeature
    finalAns$padjustGene <- ansGenes$fdrByFeature
    finalAns
}

#'
#' Mapping of chromosome names
#'
#' @param fds FraseRDataSet
#' @param style The style of the chromosome names.
#' @param ... Further parameters passed to GenomeInfoDb::mapSeqlevels().
#' 
#' @return FraseRDataSet
#' @examples
#'
#' fds <- createTestFraseRDataSet()
#' 
#' seqlevels(fds)
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


#'
#' Aberrant splicing events
#'
#' Use cutoffs to define splicing outlier events.
#'
#' @param fds FraseRDataSet object
#' @param type Splicing type
#' @param padjCutoff Adjusted p value cutoff to be used or NA if not requested
#' @param deltaPsiCutoff Delta PSI cutoff to be used or NA if not requested
#' @param zScoreCutoff Z score cutoff to be used or NA if not requested
#' @param by By default \code{none} which means no grouping. But if 
#'              \code{sample} or \code{feature} is specified the sum by 
#'              sample or feature is returned
#' @param aggregate If TRUE the returned object is based on the grouped 
#'              features
#' @param ... Further arguments can be passed to the method. If "zscores", 
#'              "padjVals" or "dPsi" is given, the values of those arguments
#'              are used to define the aberrant events.
#' @return Either a of logical values of size introns/genes x samples if "by" 
#'              is NA or a vector with the number of aberrant events per 
#'              sample or feature depending on the vaule of "by"
#' 
#' @examples
#' # get data
#' fds <- createTestFraseRDataSet()
#' 
#' # get aberrant events per sample: on the example data, nothing is aberrant
#' # based on the adjusted p-value
#' aberrant(fds, type="psi5", by="sample")
#' 
#' # use zScoreCutoff instead
#' aberrant(fds, type="psi5", by="sample", zScoreCutoff=2, padjCutoff=NA)
#' 
#' # get aberrant events per gene (first annotate gene symbols)
#' fds <- annotateRanges(fds)
#' aberrant(fds, type="psi5", by="feature", zScoreCutoff=2, padjCutoff=NA,
#'         aggregate=TRUE)
#' aberrant(fds, type="psi5", zScoreCutoff=2, padjCutoff=NA, aggregate=TRUE)
#'         
#' # find aberrant junctions/splice sites
#' aberrant(fds, type="psi5", by="none")
#' aberrant(fds, type="psi5")
#'
#' @export
#'
aberrant <- function(fds, type=currentType(fds), padjCutoff=0.05,
                    deltaPsiCutoff=0.3, zScoreCutoff=NA, minCoverage=5,
                    by=c("none", "sample", "feature"), aggregate=FALSE, ...){

    checkNaAndRange(zScoreCutoff,   min=0, max=Inf, na.ok=TRUE)
    checkNaAndRange(padjCutoff,     min=0, max=1,   na.ok=TRUE)
    checkNaAndRange(deltaPsiCutoff, min=0, max=1,   na.ok=TRUE)
    by <- match.arg(by)

    dots <- list(...)
    if("n" %in% names(dots)){
        n <- dots[['n']]
    } else {
        n <- N(fds, type=type)
    }
    if("zscores" %in% names(dots)){
        zscores <- dots[['zscores']]
    } else {
        zscores <- zScores(fds, type=type)
    }
    if("padjVals" %in% names(dots)){
        padj <- dots[['padjVals']]
    } else {
        padj <- padjVals(fds, type=type)
    }
    if("dPsi" %in% names(dots)){
        dpsi <- dots[['dPsi']]
    } else {
        dpsi <- deltaPsiValue(fds, type=type)
    } 
    
    
    # create cutoff matrix
    goodCutoff <- matrix(TRUE, nrow=nrow(zscores), ncol=ncol(zscores),
            dimnames=dimnames(zscores))
    if("hgnc_symbol" %in% colnames(mcols(fds, type=type)) &
                nrow(mcols(fds, type=type)) == nrow(goodCutoff)){
        rownames(goodCutoff) <- mcols(fds, type=type)[,"hgnc_symbol"]
    } else if(isTRUE(aggregate)){
        stop("Please provide hgnc symbols to compute gene p values!")
    }
    
    # check each cutoff if in use (not NA)
    if(!is.na(minCoverage)){
        goodCutoff <- goodCutoff & as.matrix(n >= minCoverage)
    }
    if(!is.na(zScoreCutoff)){
        goodCutoff <- goodCutoff & as.matrix(abs(zscores) > zScoreCutoff)
    }
    if(!is.na(deltaPsiCutoff)){
        goodCutoff <- goodCutoff & as.matrix(abs(dpsi) > deltaPsiCutoff)
    }
    if(!is.na(padjCutoff)){
        goodCutoff <- goodCutoff & as.matrix(padj < padjCutoff)
    }
    goodCutoff[is.na(goodCutoff)] <- FALSE
    
    # check if we should go for aggregation
    # TODO to speed it up we only use any hit within a feature
    # but should do a holm's + BY correction per gene and genome wide
    if(isTRUE(aggregate)){
        goodCutoff <- as.matrix(data.table(goodCutoff, keep.rownames=TRUE)[,
                as.data.table(t(colAnys(as.matrix(.SD)))), by=rn][,-1])
        rownames(goodCutoff) <- unique(mcols(fds, type=type)[,"hgnc_symbol"])
        colnames(goodCutoff) <- colnames(zscores)
    }
    
    # return results
    if(by == "feature"){
        return(rowSums(goodCutoff))
    }
    if(by == "sample"){
        return(colSums(goodCutoff))
    }
    return(goodCutoff)
}
