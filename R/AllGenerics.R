asSE <- function(x){
    if(is(x, "RangedSummarizedExperiment")){
        return(as(x, "RangedSummarizedExperiment"))
    }
    as(x, "SummarizedExperiment")
}

asFDS <- function(x){
    return(as(x, "FraserDataSet"))
}

#'
#' @title Getter/Setter methods for the FraserDataSet
#'
#' The following methods are getter and setter methods to extract or set
#' certain values of a FraserDataSet object. 
#' 
#' \code{samples} sets or gets the sample IDs; \code{condition} ;
#' \code{}
#' \code{nonSplicedReads} return a RangedSummarizedExperiment object 
#' containing the counts for the non spliced reads overlapping splice 
#' sites in the fds.
#' \code{}
#'
#' @param object A FraserDataSet object.
#' @param value The new value that should replace the current one.
#' @param x A FraserDataSet object.
#' @param type The psi type (psi3, psi5 or theta)
#' @return Getter method return the respective current value.
#' @examples
#' fds <- createTestFraserDataSet()
#' samples(fds)
#' samples(fds) <- 1:dim(fds)[2]
#' condition(fds)
#' condition(fds) <- 1:dim(fds)[2]
#' bamFile(fds) # file.paths or objects of class BamFile
#' bamFile(fds) <- file.path("bamfiles", samples(fds), "rna-seq.bam")
#' name(fds)
#' name(fds) <- "My Analysis"
#' workingDir(fds)
#' workingDir(fds) <- tempdir()
#' strandSpecific(fds)
#' strandSpecific(fds) <- TRUE
#' strandSpecific(fds) <- "reverse"
#' strandSpecific(fds)
#' scanBamParam(fds)
#' scanBamParam(fds) <- ScanBamParam(mapqFilter=30)
#' nonSplicedReads(fds)
#' rowRanges(fds)
#' rowRanges(fds, type="theta")
#' mcols(fds, type="psi5")
#' mcols(fds, type="theta")
#' seqlevels(fds)
#' seqlevels(mapSeqlevels(fds, style="UCSC"))
#' seqlevels(mapSeqlevels(fds, style="Ensembl"))
#' seqlevels(mapSeqlevels(fds, style="dbSNP"))
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @author Ines Scheller \email{scheller@@in.tum.de}
#'
#' @rdname fds-methods
#' @name fds-methods
NULL


#' @rdname fds-methods
#' @export
setMethod("samples", "FraserDataSet", function(object) {
    return(as.character(colData(object)[,"sampleID"]))
})

#' @rdname fds-methods
#' @export
setReplaceMethod("samples", "FraserDataSet", function(object, value) {
    colData(object)[,"sampleID"] <- as.character(value)
    rownames(colData(object)) <- colData(object)[,"sampleID"]
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("condition", "FraserDataSet", function(object) {
    if("condition" %in% colnames(colData(object))){
        return(colData(object)[,"condition"])
    }
    return(seq_len(dim(object)[2]))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("condition", "FraserDataSet", function(object, value) {
    colData(object)[,"condition"] <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("bamFile", "FraserDataSet", function(object) {
    bamFile <- colData(object)[,"bamFile"]
    if(all(vapply(bamFile, is, class2="BamFile", FUN.VALUE=logical(1)))){
        bamFile <- vapply(bamFile, path, "")
    }
    return(bamFile)
})

#' @export
#' @rdname fds-methods
setReplaceMethod("bamFile", "FraserDataSet", function(object, value) {
    colData(object)[,"bamFile"] <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("name", "FraserDataSet", function(object) {
    return(slot(object, "name"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("name", "FraserDataSet", function(object, value) {
    slot(object, "name") <- value
    validObject(object)
    return(object)
})

#' @export
#' @rdname fds-methods
setMethod("workingDir", "FraserDataSet", function(object) {
    return(slot(object, "workingDir"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("workingDir", "FraserDataSet", function(object, value) {
    slot(object, "workingDir") <- value
    validObject(object)
    return(object)
})

#' @export
#' @rdname fds-methods
setMethod("strandSpecific", "FraserDataSet", function(object) {
    return(slot(object, "strandSpecific"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("strandSpecific", "FraserDataSet", function(object, value) {
    if(is.logical(value)){
        value <- as.integer(value)
    }
    if(is.character(value)){
        value <- switch(tolower(value),
                        'no' = 0L,
                        'unstranded' = 0L,
                        'yes' = 1L,
                        'stranded' = 1L,
                        'reverse' = 2L,
                        -1L)
    }
    slot(object, "strandSpecific") <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("pairedEnd", "FraserDataSet", function(object) {
    if(!("pairedEnd" %in% colnames(colData(object)))){
        return(rep(FALSE, length(samples(object))))
    }
    pairedEnd <- colData(object)[,"pairedEnd"]
    return(pairedEnd)
})

#' @export
#' @rdname fds-methods
setReplaceMethod("pairedEnd", "FraserDataSet", function(object, value) {
    if(any(!is.logical(value))){
        warning("Value need to be logical. Converting to logical.")
    }
    colData(object)[,"pairedEnd"] <- as.logical(value)
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("scanBamParam", "FraserDataSet", function(object) {
    return(slot(object, "bamParam"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("scanBamParam", "FraserDataSet", function(object, value) {
    slot(object, "bamParam") <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("nonSplicedReads", "FraserDataSet", function(object){
    return(slot(object, "nonSplicedReads"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("nonSplicedReads", "FraserDataSet", function(object, value){
    slot(object, "nonSplicedReads") <- value
    validObject(object)
    return(object)
})

#'
#' Subsetting by indices for junctions
#'
#' Providing subsetting by indices through the single-bracket operator
#'
#' @param x A \code{FraserDataSet} object
#' @param i A integer vector to subset the rows/ranges
#' @param j A integer vector to subset the columns/samples
#' @param by a character (j or ss) definig if we subset by
#'             junctions or splice sites
#' @return A subsetted \code{FraserDataSet} object
#' @examples
#'     fds <- createTestFraserDataSet()
#'     fds[1:10,2:3]
#'     fds[,samples(fds) %in% c("sample1", "sample2")]
#'     fds[1:10,by="ss"]
#'
#' @rdname subset
subset.FRASER <- function(x, i, j, by=c("j", "ss")){
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
    subX <- as(as(x, "RangedSummarizedExperiment")[i,j], "FraserDataSet")

    # create new FraserDataSet object
    newx <- new("FraserDataSet",
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
setMethod("[", c("FraserDataSet", "ANY", "ANY"), subset.FRASER)


#'
#' Returns the assayNames of FRASER
#' 
#' @param x FraserDataSet
#' 
#' @return Character vector
#' @export
setMethod("assayNames", "FraserDataSet", function(x) {
    return(c(
        assayNames(asSE(x)),
        assayNames(nonSplicedReads(x))
    ))
})


FraserDataSet.assays.replace_pre <-
            function(x, withDimnames=TRUE, HDF5=TRUE, type=NULL, ..., value){
    if(any(names(value) == "")) stop("Name of an assay can not be empty!")
    if(any(duplicated(names(value)))) stop("Assay names need to be unique!")
    if(is.null(type)){
        type <- names(value)
    }
    # make sure all slots are HDF5
    if(isTRUE(HDF5)){
        for(i in seq_along(value)){
            if(!any(class(value[[i]]) %in% c("HDF5Matrix", "DelayedMatrix")) ||
                tryCatch(!is.character(path(value[[i]])), 
                            error=function(e){TRUE})){
                
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

    # assign new assays only for non split
    # and do the rest in the specific functions
    assays(nonSplicedReads(x), withDimnames=withDimnames, ...) <- sslots
    
    return(list(x=x, value=jslots))
}

FraserDataSet.assays.replace_r40 <- function(x, withDimnames=TRUE, HDF5=TRUE, 
                type=NULL, ..., value){
    ans <- FraserDataSet.assays.replace_pre(x, withDimnames=withDimnames, 
            HDF5=HDF5, type=type, ..., value=value)
    
    # retrieve adapted objects and set final assays on main SE object
    value <- ans[['value']]
    x <- ans[['x']]
    x <- callNextMethod()
    
    # validate and return
    validObject(x)
    return(x)
}

FraserDataSet.assays.replace_r36 <- function(x, ..., HDF5=TRUE, type=NULL, 
                withDimnames=TRUE, value){
    ans <- FraserDataSet.assays.replace_pre(x, withDimnames=withDimnames, 
            HDF5=HDF5, type=type, ..., value=value)
    
    # retrieve adapted objects and set final assays on main SE object
    value <- ans[['value']]
    x <- ans[['x']]
    x <- callNextMethod()
    
    # validate and return
    validObject(x)
    return(x)
}

FraserDataSet.assays.set_r40 <- function(x, withDimnames=TRUE, ...){
    return(c(
        assays(asSE(x), withDimnames=withDimnames, ...),
        assays(nonSplicedReads(x), withDimnames=withDimnames, ...)
    ))
}

FraserDataSet.assays.set_r36 <- function(x, ..., withDimnames=TRUE){
    FraserDataSet.assays.set_r40(x=x, ..., withDimnames=withDimnames)
}

if(compareVersion(package.version("SummarizedExperiment"), "1.17-2") < 0){
    FraserDataSet.assays.set <- FraserDataSet.assays.set_r36
    FraserDataSet.assays.replace <- FraserDataSet.assays.replace_r36
} else {
    FraserDataSet.assays.set <- FraserDataSet.assays.set_r40
    FraserDataSet.assays.replace <- FraserDataSet.assays.replace_r40
}

#'
#' Returns the assay for the given name/index of the FraserDataSet
#'
#' @param x FraserDataSet
#' @param ... Parameters passed on to SummarizedExperiment::assays() 
#' @param withDimnames Passed on to SummarizedExperiment::assays() 
#' @param HDF5 Logical value indicating whether the assay should be stored as 
#' a HDF5 file. 
#' @param type The psi type. 
#' @param value The new value to which the assay should be set. 
#'
#' @return (Delayed) matrix.
#' @export
setMethod("assays", "FraserDataSet", FraserDataSet.assays.set)

#' @rdname assays-FraserDataSet-method
#' @export
setReplaceMethod("assays", c("FraserDataSet", "SimpleList"),
        FraserDataSet.assays.replace)

#' @rdname assays-FraserDataSet-method
#' @export
setReplaceMethod("assays", c("FraserDataSet", "list"),
        FraserDataSet.assays.replace)

#' @rdname assays-FraserDataSet-method
#' @export
setReplaceMethod("assays", c("FraserDataSet", "DelayedMatrix"),
        FraserDataSet.assays.replace)


#'
#' retrieve the length of the object (aka number of junctions)
#' 
#' @param x FraserDataSet
#'
#' @return Length of the object.
#' @export
setMethod("length", "FraserDataSet", function(x) callNextMethod())

#' @rdname fds-methods 
#' @export
FRASER.mcols.get <- function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j"){
        return(mcols(asSE(x), ...))
    }
    mcols(nonSplicedReads(x), ...)
}
FRASER.mcols.replace <- function(x, type=NULL, ..., value){
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
setMethod("mcols", "FraserDataSet", FRASER.mcols.get)
setReplaceMethod("mcols", "FraserDataSet", FRASER.mcols.replace)

#' @rdname fds-methods
#' @export
FRASER.rowRanges.get <- function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j")  return(callNextMethod())
    if(type=="ss") return(rowRanges(nonSplicedReads(x), ...))
}
FRASER.rowRanges.replace <- function(x, type=NULL, ..., value){
    type <- checkReadType(x, type)
    if(type=="j") return(callNextMethod())
    rowRanges(nonSplicedReads(x), ...) <- value
    return(x)
}

setMethod("rowRanges", "FraserDataSet", FRASER.rowRanges.get)
setReplaceMethod("rowRanges", "FraserDataSet", FRASER.rowRanges.replace)

#'
#' Getter/setter for count data
#'
#' @param fds,object FraserDataSet
#' @param type The psi type.
#' @param side "ofInterest" for junction counts, "other" for sum of counts of 
#' all other junctions at the same donor site (psi5) or acceptor site (psi3), 
#' respectively. 
#' @param value An integer matrix containing the counts.
#' @param ... Further parameters that are passed to assays(object,...)
#' 
#' @return FraserDataSet
#' @rdname counts
#' @examples 
#'  fds <- createTestFraserDataSet()
#'  
#'  counts(fds, type="psi5", side="ofInterest")
#'  counts(fds, type="psi5", side="other")
#'  head(K(fds, type="psi3"))
#'  head(N(fds, type="theta"))
#'  
setMethod("counts", "FraserDataSet", function(object, type=NULL,
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
    type <- whichPSIType(type)
    if(length(type) == 0 | length(type) > 1){
        stop(paste0("Please provide a correct psi type: psi5, psi3, or ",
                    "theta. Not the given one: '", type, "'."))
    }
    aname <- paste0("rawOtherCounts_", type)
    if(!aname %in% assayNames(object)){
        stop(paste0("Missing rawOtherCounts for type '", type, "'.",
                " Please calculate PSIValues first. And then try again."))
    }
    return(assays(object)[[aname]])
})

#'
#' setter for count data
#' 
#' @rdname counts
setReplaceMethod("counts", "FraserDataSet", function(object, type=NULL,
                    side=c("ofInterest", "otherSide"), ..., value){
    side <- match.arg(side)

    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
    } else {
        type <- whichPSIType(type)
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
#' Mapping of chromosome names
#'
#' @param fds FraserDataSet
#' @param style The style of the chromosome names.
#' @param ... Further parameters. For mapSeqLevels: further parameters 
#'     passed to GenomeInfoDb::mapSeqlevels().
#' 
#' @rdname fds-methods
#' @export
mapSeqlevels <- function(fds, style="UCSC", ...){
    
    mappings <- na.omit(GenomeInfoDb::mapSeqlevels(seqlevels(fds), style, ...))
    # fix missing names() when fds has only a single chromosome
    if(is.null(names(mappings))){ 
        names(mappings) <- seqlevels(fds)
    }
    
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
#' retrieve a single sample result object
#' @noRd
resultsSingleSample <- function(sampleID, gr, pvals, padjs, zscores, 
                                psivals, rawCts, rawTotalCts, deltaPsiVals, 
                                psiType, rowMeansK, rowMeansN, aberrant, 
                                aggregate, rho, pvalsGene=NULL, padjsGene=NULL, 
                                additionalColumns, filters=list()){
    
    # if gene level results, find the most aberrant junction per gene first
    if(isTRUE(aggregate)){
        goodGenes <- rownames(aberrant)[aberrant[,sampleID]]
        geneJunctions <- findJunctionsForAberrantGenes(gr=gr, genes=goodGenes, 
                            pvals=pvals[,sampleID], 
                            dpsi=deltaPsiVals[,sampleID],
                            minCount=rawTotalCts[,sampleID], rho=rho,
                            filters=filters)    
        goodCut <- rep(FALSE, nrow(pvals))
        goodCut[geneJunctions] <- TRUE
    } else{
        goodCut <- aberrant[,sampleID]
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
    if("hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(ans)$hgncSymbol <- Rle(mcols(gr[goodCut])$hgnc_symbol)
    }
    if("other_hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(ans)$addHgncSymbols <- Rle(mcols(gr[goodCut])$other_hgnc_symbol)
    }
    mcols(ans)$type            <- Rle(psiType)
    mcols(ans)$pValue          <- signif(pvals[goodCut,sampleID], 5)
    mcols(ans)$padjust         <- signif(padjs[goodCut,sampleID], 5)
    mcols(ans)$zScore          <- Rle(round(zscores[goodCut,sampleID], 2))
    mcols(ans)$psiValue        <- Rle(round(psivals[goodCut,sampleID], 2))
    mcols(ans)$deltaPsi        <- Rle(round(deltaPsiVals[goodCut,sampleID], 2))
    mcols(ans)$meanCounts      <- Rle(round(rowMeansK[goodCut], 2))
    mcols(ans)$meanTotalCounts <- Rle(round(rowMeansN[goodCut], 2))
    mcols(ans)$counts          <- Rle(rawCts[goodCut, sampleID])
    mcols(ans)$totalCounts     <- Rle(rawTotalCts[goodCut, sampleID])
    
    if(isTRUE(aggregate)){
        mcols(ans)$pValueGene  <- signif(pvalsGene[goodCut,sampleID], 5)
        mcols(ans)$padjustGene <- signif(padjsGene[goodCut,sampleID], 5)
    }
    
    if(!is.null(additionalColumns)){
        for(column in additionalColumns){
            mcols(ans)[,column] <- Rle(mcols(gr[goodCut])[,column])
        }
    }
    
    return(ans[order(mcols(ans)$pValue)])
}

FRASER.results <- function(object, sampleIDs, fdrCutoff, zscoreCutoff, 
                            dPsiCutoff, minCount, rhoCutoff, psiType, 
                            maxCols=20, aggregate=FALSE,
                            BPPARAM=bpparam(), additionalColumns=NULL){
    
    stopifnot(is(object, "FraserDataSet"))
    stopifnot(all(sampleIDs %in% samples(object)))
    
    resultsls <- bplapply(psiType, BPPARAM=BPPARAM, function(type){
        message(date(), ": Collecting results for: ", type)
        currentType(object) <- type
        gr <- rowRanges(object, type=type)
        
        # first get row means
        rowMeansK <- rowMeans(K(object, type=type))
        rowMeansN <- rowMeans(N(object, type=type))
        
        # then iterate by chunk
        chunkCols <- getMaxChunks2Read(fds=object, assayName=type, max=maxCols)
        sampleChunks <- getSamplesByChunk(fds=object, sampleIDs=sampleIDs,
                                            chunkSize=chunkCols)
        
        ans <- lapply(seq_along(sampleChunks), function(idx){
            message(date(), ": Process chunk: ", idx, " for: ", type)
            sc <- sampleChunks[[idx]]
            tmp_x <- object[,sc]
            
            # extract values
            rawCts       <- as.matrix(K(tmp_x))
            rawTotalCts  <- as.matrix(N(tmp_x))
            pvals        <- as.matrix(pVals(tmp_x))
            padjs        <- as.matrix(padjVals(tmp_x))
            zscores      <- as.matrix(zScores(tmp_x))
            psivals      <- as.matrix(assay(tmp_x, type))
            muPsi        <- as.matrix(predictedMeans(tmp_x))
            psivals_pc   <- (rawCts + pseudocount()) /
                (rawTotalCts + 2*pseudocount())
            deltaPsiVals <- psivals_pc - muPsi
            rho          <- rho(tmp_x, type)
            aberrant     <- aberrant.FRASER(tmp_x, type=type, 
                                                padjCutoff=fdrCutoff, 
                                                zScoreCutoff=zscoreCutoff, 
                                                deltaPsiCutoff=dPsiCutoff, 
                                                minCount=minCount, 
                                                rhoCutoff=rhoCutoff,
                                                aggregate=aggregate)
            if(isTRUE(aggregate)){
                pvalsGene    <- as.matrix(pVals(tmp_x, level="gene"))
                padjsGene    <- as.matrix(padjVals(tmp_x, level="gene"))
            } else{
                pvalsGene    <- NULL
                padjsGene    <- NULL
            }
            
            if(length(sc) == 1){
                colnames(pvals) <- sc
                colnames(padjs) <- sc
                colnames(zscores) <- sc
                colnames(deltaPsiVals) <- sc
            }
            
            # create result table
            sampleRes <- lapply(sc,
                                resultsSingleSample, gr=gr, pvals=pvals, 
                                padjs=padjs, zscores=zscores, psiType=type, 
                                psivals=psivals, deltaPsiVals=deltaPsiVals, 
                                rawCts=rawCts, rawTotalCts=rawTotalCts, 
                                rowMeansK=rowMeansK, rowMeansN=rowMeansN, 
                                aberrant=aberrant, aggregate=aggregate,
                                rho=rho, 
                                pvalsGene=pvalsGene, padjsGene=padjsGene,
                                additionalColumns=additionalColumns,
                                filters=list(dpsi=dPsiCutoff, 
                                            minCount=minCount, 
                                            rho=rhoCutoff))
            
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
#' Extracting results and aberrant splicing events
#'
#' The result function extracts the results from the given analysis object
#' based on the given options and cutoffs. The aberrant function extracts 
#' aberrant splicing events based on the given cutoffs.
#'
#' @param object A \code{\link{FraserDataSet}} object
#' @param sampleIDs A vector of sample IDs for which results should be 
#' retrieved
#' @param padjCutoff The FDR cutoff to be applied or NA if not requested.
#' @param zScoreCutoff The z-score cutoff to be applied or NA if not requested.
#' @param deltaPsiCutoff The cutoff on delta psi or NA if not requested.
#' @param minCount The minimum count value of the total coverage of an intron 
#' to be considered as significant.
#' result
#' @param psiType The psi types for which the results should be retrieved.
#' @param additionalColumns Character vector containing the names of additional 
#' columns from mcols(fds) that should appear in the result table 
#' (e.g. ensembl_gene_id). Default is \code{NULL}, so no additional columns 
#' are included. 
#' @param BPPARAM The BiocParallel parameter.
#' @param res Result as created with \code{results()}
#' @param geneColumn The name of the column in \code{mcols(res)} that contains 
#'     the gene symbols.   
#' @param method The p.adjust method that is being used to adjust p values per
#'     sample.
#' @param type Splicing type (psi5, psi3 or theta)
#' @param by By default \code{none} which means no grouping. But if 
#'              \code{sample} or \code{feature} is specified the sum by 
#'              sample or feature is returned
#' @param aggregate If TRUE the returned object is based on the grouped 
#'              features
#' @param ... Further arguments can be passed to the method. If "zscores", 
#'              "padjVals" or "dPsi" is given, the values of those arguments
#'              are used to define the aberrant events.
#'
#' @return For \code{results}: GRanges object containing significant results.
#'     For \code{aberrant}: Either a of logical values of size 
#'     introns/genes x samples if "by" is NA or a vector with the 
#'     number of aberrant events per sample or feature depending on 
#'     the vaule of "by"
#' 
#' @rdname results
#' @examples
#' # get data, fit and compute p-values and z-scores
#' fds <- createTestFraserDataSet()
#' 
#' # extract results: for this example dataset, z score cutoff of 2 is used to
#' # get at least one result and show the output
#' res <- results(fds, padjCutoff=NA, zScoreCutoff=3, deltaPsiCutoff=0.05)
#' res
#' 
#' # aggregate the results by genes (gene symbols need to be annotated first 
#' # using annotateRanges() function)
#' resultsByGenes(res)
#'
#' # get aberrant events per sample: on the example data, nothing is aberrant
#' # based on the adjusted p-value
#' aberrant(fds, type="psi5", by="sample")
#' 
#' # get aberrant events per gene (first annotate gene symbols)
#' fds <- annotateRangesWithTxDb(fds)
#' aberrant(fds, type="psi5", by="feature", zScoreCutoff=2, padjCutoff=NA,
#'         aggregate=TRUE)
#'         
#' # find aberrant junctions/splice sites
#' aberrant(fds, type="psi5")
#' @export
setMethod("results", "FraserDataSet", function(object, 
                    sampleIDs=samples(object), padjCutoff=0.05,
                    zScoreCutoff=NA, deltaPsiCutoff=0.3,
                    rhoCutoff=0.1, aggregate=FALSE,
                    minCount=5, psiType=c("psi3", "psi5", "theta"),
                    additionalColumns=NULL, BPPARAM=bpparam(), ...){
    FRASER.results(object=object, sampleIDs=sampleIDs, fdrCutoff=padjCutoff,
            zscoreCutoff=zScoreCutoff, dPsiCutoff=deltaPsiCutoff, 
            rhoCutoff=rhoCutoff, minCount=minCount, 
            psiType=match.arg(psiType, several.ok=TRUE),
            aggregate=aggregate,
            additionalColumns=additionalColumns, BPPARAM=BPPARAM)
})

#' @rdname results
#' @export
resultsByGenes <- function(res, geneColumn="hgncSymbol"){
    # sort by gene pvalue
    res <- res[order(res$pValueGene)]

    # extract subset
    if(is(res, "GRanges")){
        ans <- as.data.table(mcols(res)[,c(geneColumn, "sampleID")])
        colnames(ans) <- c("features", "sampleID")
    } else {
        ans <- featureNames <- res[,.(
                features=get(geneColumn), sampleID=sampleID)]
    }

    # remove NAs and 
    # keep only one row per gene with lowest pvalue over all psiTypes
    naIdx <- ans[,is.na(features)]
    ansNoNA <- ans[!is.na(features)]
    dupIdx <- duplicated(ansNoNA[,.(features,sampleID)])

    # get final result table
    finalAns <- res[!naIdx][!dupIdx]
    finalAns
}

aberrant.FRASER <- function(object, type=currentType(object), 
                                padjCutoff=0.05, deltaPsiCutoff=0.3, 
                                zScoreCutoff=NA, minCount=5, rhoCutoff=0.1,
                                by=c("none", "sample", "feature"), 
                                aggregate=FALSE, ...){
    
    checkNaAndRange(padjCutoff,     min=0, max=1,   scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(zScoreCutoff,   min=0, max=Inf, scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(deltaPsiCutoff, min=0, max=1,   scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(rhoCutoff,      min=0, max=1,   scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(minCount,       min=0, max=Inf, scalar=TRUE,   na.ok=TRUE)
    by <- match.arg(by)
    
    dots <- list(...)
    if("n" %in% names(dots)){
        n <- dots[['n']]
    } else {
        n <- N(object, type=type)
    }
    if("zscores" %in% names(dots)){
        zscores <- dots[['zscores']]
    } else {
        zscores <- zScores(object, type=type)
    }
    if("padjVals" %in% names(dots)){
        padj <- dots[['padjVals']]
    } else {
        # check if padj values are available for the given filters
        pvalsAvailable <- checkPadjAvailableForFilters(object, type=type,
                                                filters=list(rho=rhoCutoff), 
                                                aggregate=aggregate)
        if(isFALSE(pvalsAvailable)){
            stop("For the given filters, pvalues are not yet computed. \n", 
                    "Please compute them first by running the ",
                    "calculatePadjValues function with the requested filters.")
        }
        pvalLevel <- ifelse(isTRUE(aggregate), "gene", "site")
        padj <- padjVals(object, type=type, level=pvalLevel, 
                        filters=list(rho=rhoCutoff))
    }
    if("dPsi" %in% names(dots)){
        dpsi <- dots[['dPsi']]
    } else {
        dpsi <- deltaPsiValue(object, type=type)
    } 
    if("rhoVals" %in% names(dots)){
        rho <- dots[['rhoVals']]
    } else {
        rho <- matrix(rho(object, type=type), 
                        nrow=nrow(dpsi), ncol=ncol(dpsi))
    } 
    
    if(is.na(padjCutoff)){
        padjCutoff <- 1
    }
    
    if(isFALSE(aggregate)){
        aberrantEvents <- padj <= padjCutoff
        
        # check each cutoff if in use (not NA)
        if(!is.na(minCount)){
            aberrantEvents <- aberrantEvents & as.matrix(n >= minCount)
        }
        if(!is.na(zScoreCutoff)){
            aberrantEvents <- aberrantEvents & 
                                as.matrix(abs(zscores) > zScoreCutoff)
        }
        if(!is.na(deltaPsiCutoff)){
            aberrantEvents <- aberrantEvents & 
                                as.matrix(abs(dpsi) > deltaPsiCutoff)
        }
        if(!is.na(rhoCutoff)){
            aberrantEvents <- aberrantEvents & as.matrix(rho < rhoCutoff)
        }
        aberrantEvents[is.na(aberrantEvents)] <- FALSE
        
    } else{
        dt <- data.table(geneID=mcols(fds, type=type)$hgnc_symbol, padj)
        dt <- dt[!duplicated(dt, by="geneID") & !is.na(geneID)]
        padj <- as.matrix(dt[, -1])
        rownames(padj) <- dt[,geneID]
        
        aberrantEvents <- padj <= padjCutoff
    }
    
    return(switch(match.arg(by),
                    none = aberrantEvents,
                    sample = colSums(aberrantEvents, na.rm=TRUE),
                    feature = rowSums(aberrantEvents, na.rm=TRUE)
    ))
}

#' @rdname results
#' @export
setMethod("aberrant", "FraserDataSet", aberrant.FRASER)


