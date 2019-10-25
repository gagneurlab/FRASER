##
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for reading in data
## especially aligned RNA sequencing data
##

#'
#' @title Count RNA-seq data
#' @description This method extracts and counts the split reads and
#'             non spliced reads from a RNA bam file
#' @details TODO
#'
#' @param fds A FraseRDataSet object with all the information
#'             how and what to count
#' @param NcpuPerSample A BiocParallel param object or a positive integer
#'             to configure the parallel backend
#'             of the internal loop per sample
#' @param junctionMap A object or file containing a map of
#'             all junctions of interest across all samples
#' @param minAnchor Minimum overlap around the Donor/Acceptor for
#'             non spliced reads. Default to 5
#' @param recount if TRUE the cache is ignored and the bam file is recounted.
#'
#' @return FraseRDataSet
#' @export
#' @examples
#'   countRNAData(createTestFraseRSettings())
#'   countRNAData(createTestFraseRSettings(), 5)
countRNAData <- function(fds, NcpuPerSample=1, junctionMap=NULL, minAnchor=5,
                    recount=FALSE, BPPARAM=parallel(fds), genome=NULL){

    # Check input TODO
    stopifnot(class(fds) == "FraseRDataSet")
    stopifnot(is.numeric(NcpuPerSample) && NcpuPerSample > 0)
    stopifnot(is.numeric(minAnchor) & minAnchor >= 1)
    minAnchor <- as.integer(minAnchor)

    # load needed genomes if provided
    if(!(is.null(genome) | any(is.na(unique(genome))))){
        for(i in unique(genome)){
            library(i, character.only=TRUE)
        }
    }

    if(!is.integer(NcpuPerSample)){
        NcpuPerSample <- as.integer(NcpuPerSample)
    }
    if(!is.null(junctionMap)){
        junctionMap <- readJunctionMap(junctionMap)
    }

    # count splitreads first
    message(date(), ": Start counting the split reads ...")
    countList <- bplapply(samples(fds),
            FUN=countSplitReads,
            settings=fds,
            BPPARAM=BPPARAM,
            NcpuPerSample=NcpuPerSample,
            genome=genome,
            recount=recount)
    names(countList) <- samples(fds)
    BPPARAM_old <- setMaxThreads(BPPARAM, 10, 10)
    counts <- mergeCounts(countList, junctionMap=junctionMap,
            assumeEqual=FALSE, BPPARAM=BPPARAM_old$BPPARAM)
    BPPARAM <- setMaxThreads(BPPARAM_old$BPPARAM,
            BPPARAM_old$numWorkers, BPPARAM_old$numTasks, set=TRUE)

    # create summarized objects
    h5 <- saveAsHDF5(fds, "rawCountsJ", as.matrix(mcols(counts)[samples(fds)]))
    colnames(h5) <- samples(fds)
    splitCounts <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsJ=h5),
        rowRanges=counts[,!colnames(mcols(counts)) %in% samples(fds)]
    )
    rm(countList, counts, h5)
    gc()

    # splice site map
    splitCounts <- annotateSpliceSite(splitCounts)

    # count the retained reads
    message(date(), ": Start counting the non spliced reads ...")
    message(date(), ": In total ", length(splitCounts),
            " splice junctions are found."
    )

    # extract donor and acceptor sites
    spliceSiteCoords <- extractSpliceSiteCoordinates(splitCounts, fds)
    message(date(), ": In total ", length(spliceSiteCoords),
            " splice sites (acceptor/donor) will be counted ..."
    )

    # count non spliced reads
    countList <- bplapply(samples(fds),
            FUN=countNonSplicedReads,
            settings=fds,
            spliceSiteCoords=spliceSiteCoords,
            BPPARAM=BPPARAM,
            NcpuPerSample=NcpuPerSample,
            minAnchor=minAnchor,
            recount=recount)

    names(countList) <- samples(fds)
    siteCounts <- mergeCounts(countList, assumeEqual=TRUE)
    mcols(siteCounts)$type <- factor(countList[[1]]$type,
            levels = c("Acceptor", "Donor")
    )

    # create summarized objects
    h5 <- saveAsHDF5(fds, "rawCountsSS",
            as.matrix(mcols(siteCounts)[samples(fds)]))
    colnames(h5) <- samples(fds)
    nonSplicedCounts <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsSS=h5),
        rowRanges=siteCounts[,!colnames(mcols(siteCounts)) %in% samples(fds)]
    )
    rm(countList, siteCounts, h5)
    gc()

    # create final FraseR dataset
    fds <- new("FraseRDataSet",
        splitCounts,
        name            = name(fds),
        method          = method(fds),
        parallel        = parallel(fds),
        bamParam        = scanBamParam(fds),
        strandSpecific  = strandSpecific(fds),
        workingDir      = workingDir(fds),
        nonSplicedReads = nonSplicedCounts,
        metadata        = metadata(fds)
    )

    # save it so the FraseR object also gets saved
    fds <- saveFraseRDataSet(fds)
    gc()

    # return it
    return(fds)
}

#'
#' extracts the chromosomes within the given bamFile
#' @noRd
extractChromosomes <- function(bamFile){
    names(scanBamHeader(bamFile)[[bamFile]]$target)
}


#'
#' returns the name of the cache file if caching is enabled
#' for the given sample
#' @noRd
getSplitCountCacheFile <- function(sampleID, settings){

    # cache folder
    cachedir <- file.path(workingDir(settings), "cache", "splitCounts")
    if(!dir.exists(cachedir)){
        dir.create(cachedir, recursive=TRUE)
    }

    # file name
    filename <- paste0("splitCounts-", sampleID, ".RDS")

    # return it
    return(file.path(cachedir, filename))
}


#'
#' count all split reads in a bam file
#' @noRd
countSplitReads <- function(sampleID, settings, NcpuPerSample, genome, recount){
    message(date(), ": Count split reads for sample: ", sampleID)
    bamFile <- bamFile(settings[,samples(settings) == sampleID])[[1]]

    # check cache if available
    cacheFile <- getSplitCountCacheFile(sampleID, settings)
    if(isFALSE(recount) && !is.null(cacheFile) && file.exists(cacheFile)){
        cache <- readRDS(cacheFile)
        bamWhich <- unlist(bamWhich(scanBamParam(settings)))
        if(length(bamWhich) > 0){
            userTargetGR <- GRanges(seqnames=names(unlist(bamWhich)),
                    ranges=unlist(bamWhich), strand="*")
            from <- unique(from(findOverlaps(cache, userTargetGR)))
            cache <- cache[from]
        }
        if(length(cache) > 0){
            return(checkSeqLevelStyle(cache, settings, sampleID,
                    sampleSpecific=FALSE))
        }
    }

    # parallelize over chromosomes
    chromosomes <- extractChromosomes(bamFile)

    if(length(genome) > 1){
        genome <- genome[sampleID]
    }

    # extract the counts per chromosome
    countsList <- bplapply(chromosomes,
            FUN=countSplitReadsPerChromosome,
            bamFile=bamFile,
            settings=settings,
            BPPARAM=MulticoreParam(NcpuPerSample, length(chromosomes)),
            genome=genome)

    # sort and merge the results befor returning/saving
    countsGR <- sort(unlist(GRangesList(countsList)))
    saveRDS(countsGR, cacheFile)

    return(checkSeqLevelStyle(countsGR, settings, sampleID, sampleSpecific=FALSE))
}


#'
#' counting the split reads per chromosome
#' @noRd
countSplitReadsPerChromosome <- function(chromosome, bamFile, settings, genome){
    # restrict to the chromosome only
    which=GRanges(
        seqnames=chromosome,
        ranges=IRanges(0, 536870912)
    )
    param <- mergeBamParams(bamParam=scanBamParam(settings), which=which)
    if(is.null(param)){
        return(GRanges())
    }

    # get reads from bam file
    galignment <- readGAlignments(bamFile, param=param)

    # remove the strand information if unstranded data
    if(!strandSpecific(settings)){
        strand(galignment) <- "*"
    }

    # dont count if there is nothing to count
    if(length(galignment) == 0){
        return(GRanges())
    }

    # get the junction positions and their counts

    jc <- summarizeJunctions(galignment, genome=genome)

    ans <- jc[,"score"]
    colnames(mcols(ans)) <- "count"

    # set predicted strand if present or set it to + if NA
    if(!strandSpecific(settings) & !is.null(genome) & length(ans) > 0){
        strand(ans) <- jc$intron_strand
        ans$intron_motif <- jc$intron_motif

        # set remaining unknown junction to plus strand
        # (its 50/50 that we are wrong)
        strand(ans)[jc$intron_strand == "*"] <- "+"
    }

    # sort it and return the GRange object
    sort(ans)
}


#'
#' merge a ScanBamParam object with a given which object (GRange)
#' @noRd
mergeBamParams <- function(bamParam, which, override=FALSE){
    # the chromosome is alrways only one
    chromosome <- as.character(unique(seqnames(which)))

    # just take the which argument of no ranges are specified
    if(length(bamWhich(bamParam)) == 0 | override){
        bamWhich(bamParam) <- which
    } else {
        if(is.null(bamWhich(bamParam)[[chromosome]])){
            return(NULL)
        }
        # only take the ranges overlapping with the given once
        ov <- findOverlaps(bamWhich(bamParam)[[chromosome]], ranges(which))
        bamWhich(bamParam)[[chromosome]] <-
                bamWhich(bamParam)[[chromosome]][from(ov)]
    }
    return(bamParam)
}


#'
#' merge multiple samples into one SummarizedExperiment object
#' @noRd
mergeCounts <- function(countList, junctionMap=NULL, assumeEqual=FALSE,
                BPPARAM=SerialParam()){

    # prepare range object
    sample_names <- names(countList)

    if(assumeEqual){
        ranges <- sort(unique(GRanges(countList[[1]])))
        mcols(ranges)$count <- NULL
        names(ranges) <- NULL

        message(paste(date(), ": Fast merging of counts ..."))
        sample_counts <- lapply(countList, function(gr) mcols(gr)[["count"]] )
    } else {
        countList <- uniformSeqInfo(countList)
        countgr <- GRangesList(countList)
        if(!is.null(junctionMap)){
            stopifnot(class(junctionMap) == "GRanges")
            countgr <- c(countgr, GRangesList(junctionMap))
        }
        ranges <- sort(unique(unlist(countgr)))
        mcols(ranges)$count <- NULL
        names(ranges) <- NULL

        message(paste(date(), ": count ranges need to be merged ..."))
        # merge each sample counts into the combined range object
        sample_counts <- bplapply(countList, ranges = ranges, BPPARAM=BPPARAM,
            FUN = function(gr, ranges){
                suppressPackageStartupMessages(require(FraseR))

                # init with 0 since we did not find any read
                # for this sample supporting a given site
                sample_count <- integer(length(ranges))

                # get overlap and add counts to the corresponding ranges
                overlaps <- findOverlaps(gr, ranges, type = "equal")
                sample_count[overlaps@to] <- mcols(gr)$count

                return(sample_count)
            }
        )
    }

    # merge it with the type columen and add it to the range object
    counts <- DataFrame(sample_counts)
    colnames(counts) <- sample_names
    mcols(ranges) <- cbind(mcols(ranges), counts)

    # return the object
    return(ranges)
}


#'
#' returns the name of the cache file if caching is enabled for the given sample
#' @noRd
getNonSplicedCountCacheFile <- function(sampleID, settings){

    # cache folder
    cachedir <- file.path(workingDir(settings), "cache", "nonSplicedCounts")
    if(!dir.exists(cachedir)){
        dir.create(cachedir, recursive=TRUE)
    }

    # file name
    filename <- paste0("nonSplicedCounts-", sampleID, ".RDS")

    # return it
    return(file.path(cachedir, filename))
}


#'
#' creates a SAF data.table based on the given grange like object
#' @noRd
GRanges2SAF <- function(gr, minAnchor=1){
    data.table(
        GeneID  = seq_along(gr),
        Chr     = as.factor(seqnames(gr)),
        Start   = start(gr) - (minAnchor - 1),
        End     = end(gr) + (minAnchor - 1),
        Strand  = as.factor(strand(gr))
    )
}


#'
#' counts non spliced reads based on the given target (acceptor/donor) regions
#' TODO: 10k chunks hardcoded currently (needs some testing and a code to)
#' @noRd
countNonSplicedReads <- function(sampleID, spliceSiteCoords, settings,
                    NcpuPerSample=1, minAnchor, recount=FALSE){

    message(date(), ": Count non spliced reads for sample: ", sampleID)
    bamFile <- bamFile(settings[,samples(settings) == sampleID])[[1]]

    # check cache if available
    cacheFile <- getNonSplicedCountCacheFile(sampleID, settings)
    if(isFALSE(recount) & !is.null(cacheFile) && file.exists(cacheFile)){
        # check if needs to be recalculated
        cache <- try(readRDS(cacheFile), silent=TRUE)
        if(class(cache) == "try-error"){
            message(date(), ":Cache is currupt for sample:",
                    sampleID, ". Will recount it."
            )
            cache <- GRanges()
        }
        cache <- checkSeqLevelStyle(cache, settings, sampleID, FALSE)
        ov    <- findOverlaps(cache, spliceSiteCoords, type="equal")

        # we have all sites of interest cached
        if(all(1:length(spliceSiteCoords) %in% to(ov))){
            cache2return <- cache[from(ov)]
            cache2return$spliceSiteID <- spliceSiteCoords[to(ov)]$spliceSiteID
            cache2return$type <- spliceSiteCoords[to(ov)]$type
            return(cache2return)
        } else {
            message(paste("The cache file does not contain the needed",
                    "genomic positions. Adding the remaining sites to",
                    "the cache ..."
            ))
        }
    }

    # extract the counts with Rsubread
    tmp_ssc <- checkSeqLevelStyle(spliceSiteCoords, settings, sampleID, TRUE)
    anno <- GRanges2SAF(tmp_ssc, minAnchor=minAnchor)
    rsubreadCounts <- featureCounts(files=bamFile, annot.ext=anno,
            minOverlap=minAnchor*2, allowMultiOverlap=TRUE,
            checkFragLength=FALSE,
            minMQS=bamMapqFilter(scanBamParam(settings)),
            strandSpecific=as.integer(strandSpecific(settings)),

            # activating long read mode
            isLongRead=TRUE,

            # multi-mapping reads
            countMultiMappingReads=TRUE,

            # for counting only non spliced reads we skip this information
            isPairedEnd=FALSE,

            # this is important, otherwise it will sort it by name
            autosort=FALSE,
            nthreads=NcpuPerSample,
            tmpDir=file.path(workingDir(settings), "cache")
    )

    # extract results
    mcols(spliceSiteCoords)$count <- rsubreadCounts$counts[,1]
    spliceSiteCoords <- sort(spliceSiteCoords)


    # cache it if enabled
    if(!is.null(cacheFile)){
        if(exists("cache") && class(cache) == "GRanges"){
            message("Adding splice sites to cache ...")
            cache <- unique(unlist(GRangesList(spliceSiteCoords, cache)))
        } else {
            message("Saving splice site cache ...")
            cache <- spliceSiteCoords
        }
        cache <- checkSeqLevelStyle(cache, settings, sampleID, TRUE)
        saveRDS(cache, cacheFile)
    }


    # return it
    return(spliceSiteCoords)
}

#'
#' reads the given global junction map and merges it with
#' the given junction map
#'
#' @noRd
readJunctionMap <- function(junctionMap){
    if(is.null(junctionMap)){
        return(NULL)
    }

    # read global junction map
    if(isScalarCharacter(junctionMap)){
        stopifnot(file.exists(junctionMap))
        if(endsWith(junctionMap, ".RDS")){
            map <- readRDS(junctionMap)
            stopifnot(class(map) == "GRanges")
            return(map)
        }
        message(date(), ": currently not supported.")
        message(date(), ": will only use the provided junctions")
        return(NULL)
    }
    if(class(junctionMap) == "GRanges"){
        return(junctionMap)
    }
    stop("Objecttype (class) for junction map currently not supported:",
            class(junctionMap)
    )
}

#'
#' extracts the splice site coordinates from a junctions GRange object
#' @noRd
extractSpliceSiteCoordinates <- function(junctions, fds){

    if(isTRUE(strandSpecific(fds))){
        spliceSiteCoords <- unlist(GRangesList(
            extractSpliceSiteCoordsPerStrand(junctions, "+"),
            extractSpliceSiteCoordsPerStrand(junctions, "-")
        ))
    } else {
        strand(junctions) <- "*"
        spliceSiteCoords <- extractSpliceSiteCoordsPerStrand(junctions, "*")
    }

    return(unique(sort(spliceSiteCoords)))
}


#'
#' extracts the splice site coordinates per strand from a
#' given junctions GRange object
#' @noRd
extractSpliceSiteCoordsPerStrand <- function(junctions, strand){

    # get only the correct strand features
    junctions <- junctions[strand(junctions) == strand,]

    # left side (acceptor on + and donor on -)
    left <- GRanges(
            seqnames = seqnames(junctions),
            strand = strand(junctions),
            ranges = IRanges(
                    start = start(junctions) - 1,
                    end   = start(junctions)
            ),
            seqlengths = seqlengths(junctions),
            seqinfo = seqinfo(junctions),
            mcols(junctions)[c("startID")]
    )
    colnames(mcols(left)) <- "spliceSiteID"

    # right side (acceptor on - and donor on +)
    right <- GRanges(
            seqnames = seqnames(junctions),
            strand = strand(junctions),
            ranges = IRanges(
                    start = end(junctions),
                    end   = end(junctions) + 1
            ),
            seqlengths = seqlengths(junctions),
            seqinfo = seqinfo(junctions),
            mcols(junctions)[c("endID")]
    )
    colnames(mcols(right)) <- "spliceSiteID"

    # annotate donor and acceptor sites
    if(strand %in% c("+", "*")){
        mcols(left)$type  = "Donor"
        mcols(right)$type = "Acceptor"
    } else {
        mcols(left)$type  = "Acceptor"
        mcols(right)$type = "Donor"
    }

    return(unique(sort(unlist(GRangesList(left, right)))))
}

#'
#' annotates the given GRange object with unique identifier for
#' splice sites (acceptors and donors)
#' @noRd
annotateSpliceSite <- function(gr){
    message(date(), ": Create splice site indices ...")

    # convert to data.table for better handling
    dt <- GRanges2SAF(gr)

    # extract donor/acceptor annotation
    startSideDT <- dt[,.(End=Start, type="start"),by="Chr,Start,Strand"]
    endSideDT   <- dt[,.(Start=End, type="end"  ),by="Chr,End,Strand"]

    # annotate and enumerate donor/acceptor
    annotadedDT <- rbind(startSideDT, endSideDT)
    annotadedDT[,id:=1:nrow(annotadedDT)]

    # convert back to granges
    annogr <- makeGRangesFromDataFrame(annotadedDT, keep.extra.columns=TRUE)

    ids <- mclapply(c("start", "end"), mc.cores=2, function(type){
        # reduce annogr to only the specific type to prevent overlap
        annogrtmp <- annogr[annogr$type == type]
        ov <- findOverlaps(gr, annogrtmp, type=type)
        mcols(annogrtmp[to(ov)])[["id"]]
    })

    names(ids) <- paste0(c("start", "end"), "ID")
    mcols(gr)[names(ids)] <- DataFrame(ids)

    return(gr)
}

setMaxThreads <- function(BPPARAM, maxWorkers=bpworkers(BPPARAM), maxTasks=bptasks(BPPARAM), set=FALSE){
    numWorkers <- bpworkers(BPPARAM)
    numTasks <- bptasks(BPPARAM)
    try({
        if(isFALSE(set)){
            maxWorkers <- min(maxWorkers, numWorkers)
            maxTasks <- min(maxTasks, numTasks)
        }
        bpworkers(BPPARAM) <- maxWorkers
        bptasks(BPPARAM) <- maxTasks
    }, silent=TRUE)
    if(isTRUE(set)){
        return(BPPARAM)
    }
    return(list(BPPARAM=BPPARAM, numWorkers=numWorkers, numTasks=numTasks))
}
