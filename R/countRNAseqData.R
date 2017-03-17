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
#' @param settings A FraseRSetting object with all the information
#'             how and what to count
#' @param internBPPARAM A BiocParallel param object to configure the
#'             parallel backend of the internal loop
#'
#' @return FraseRDataSet
#' @export
#' @examples
#'   counRNAData(createTestFraseRSettings())
countRNAData <- function(settings, internBPPARAM=SerialParam()){

    # Check input TODO
    stopifnot(class(settings) == "FraseRSettings")
    stopifnot(is(internBPPARAM, "BiocParallelParam"))

    # count splitreads first
    message(date(), ": Start counting the split reads ...")
    countList <- bplapply(samples(settings),
            FUN=.countSplitReads,
            settings=settings,
            BPPARAM=parallel(settings),
            internBPPARAM=internBPPARAM
    )
    names(countList) <- samples(settings)
    counts <- .mergeCounts(countList, assumeEqual=FALSE, parallel(settings))

    # count the retained reads
    message(date(), ": Start counting the non spliced reads ...")
    message(date(), ": In total ", length(granges(counts)),
            " splice junctions are found."
    )

    # extract donor and acceptor sites
    spliceSiteCoords <- .extractSpliceSiteCoordinates(granges(counts), settings)
    message(date(), ": In total ", length(spliceSiteCoords),
            " splice sites (acceptor/donor) will be counted ..."
    )

    #
    # TODO use spliceSites instead of targets to reduce computation
    #
    countList <- bplapply(samples(settings),
            FUN=.countNonSplicedReads,
            settings=settings,
            spliceSiteCoords=spliceSiteCoords,
            BPPARAM=parallel(settings),
            internBPPARAM=internBPPARAM
    )
    names(countList) <- samples(settings)
    site_counts <- .mergeCounts(countList, assumeEqual=TRUE, parallel(settings))
    mcols(site_counts)$type=factor(countList[[1]]$type,
            levels = c("Acceptor", "Donor")
    )

    # create summarized objects
    splitCounts <- SummarizedExperiment(
            assays=list(rawCounts=mcols(counts)[samples(settings)]),
            rowRanges=granges(counts)
    )
    nonSplicedCounts <- SummarizedExperiment(
            assays=list(rawCounts=mcols(site_counts)[samples(settings)]),
            rowRanges=site_counts[,"type"]
    )

    # create final FraseR dataset
    fds <- FraseRDataSet(
        settings=settings,
        splitReads=splitCounts,
        nonSplicedReads=nonSplicedCounts
    )

    # return it
    return(fds)
}

##
## extracts the chromosomes within the given bamFile
##
.extractChromosomes <- function(bamFile){
    names(scanBamHeader(path(bamFile))[[path(bamFile)]]$target)
}


##
## returns the name of the cache file if caching is enabled
## for the given sample
##
.getSplitCountCacheFile <- function(sampleID, settings){
    # check if caching is enabled
    if(is.null(outputFolder(settings)) || outputFolder(settings) == "")
        return(NULL)

    # cache folder
    cachedir <- file.path(outputFolder(settings), "cache", "splitCounts")
    if(!dir.exists(cachedir)){
        dir.create(cachedir, recursive=TRUE)
    }

    # file name
    filename <- paste0("splitCounts-", sampleID, ".RDS")

    # return it
    return(file.path(cachedir, filename))
}


##
## count all split reads in a bam file
##
.countSplitReads <- function(sampleID, settings, internBPPARAM){
    suppressPackageStartupMessages(library(FraseR))
    bamFile <- bamFiles(settings[samples(settings) == sampleID])[[1]]

    # check cache if available
    cacheFile <- .getSplitCountCacheFile(sampleID, settings)
    if(!is.null(cacheFile) && file.exists(cacheFile)){
        return(readRDS(cacheFile))
    }

    # parallelize over chromosomes
    chromosomes <- .extractChromosomes(bamFile)

    # extract the counts per chromosome
    countsList <- bplapply(chromosomes,
            FUN=.countSplitReadsPerChromosome,
            bamFile=bamFile,
            settings=settings,
            BPPARAM=internBPPARAM
    )

    # sort and merge the results befor returning/saving
    countsGR <- sort(unlist(GRangesList(countsList)))
    if(!is.null(cacheFile)){
        saveRDS(countsGR, cacheFile)
    }
    return(countsGR)
}


##
## counting the split reads per chromosome
##
.countSplitReadsPerChromosome <- function(chromosome, bamFile, settings){
    # restrict to the chromosome only
    which=GRanges(
        seqnames=chromosome,
        ranges=IRanges(0, 536870912)
    )
    param <- .mergeBamParams(bamParam=settings@bamParams, which=which)
    if(is.null(param)){
        return(GRanges())
    }

    # get reads from bam file
    galignment <- readGAlignments(bamFile, param=param)

    # remove the strand information if unstranded data
    if(!settings@strandSpecific){
        strand(galignment) <- "*"
    }

    # dont count if there is nothing to count
    if(length(galignment) == 0){
        return(GRanges())
    }

    # get the junction positions and their counts
    junctions <- unlist(junctions(galignment))
    junctionsCounts <- unique(junctions)

    # dont count anything if there is nothing to count
    if(length(junctionsCounts) == 0){
        return(junctionsCounts)
    }

    # count the data
    mcols(junctionsCounts)$count <- countOverlaps(
            junctionsCounts,
            junctions,
            type = 'equal'
    )

    # sort it and return the GRange object
    return(sort(junctionsCounts))
}


##
## merge a ScanBamParam object with a given which object (GRange)
##
.mergeBamParams <- function(bamParam, which, override=FALSE){
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


##
## merge multiple samples into one SummarizedExperiment object
##
## TODO: how to init non found junctions/splice sites (0L or NA)?
##
.mergeCounts <- function(countList, assumeEqual=FALSE, BPPARAM=SerialParam()){

    # prepare range object
    sample_names <- names(countList)
    countgr <- GRangesList(countList)
    ranges <- sort(unique(unlist(countgr)))
    mcols(ranges)$count <- NULL
    names(ranges) <- NULL

    if(assumeEqual){
        message(paste(date(), ": Fast merging of counts ..."))
        sample_counts <- lapply(1:length(countgr), function(i){
                mcols(countgr[[i]])$count
        })
    } else {
        message(paste(date(), ": count ranges need to be merged ..."))
        # merge each sample counts into the combined range object
        sample_counts <- bplapply(1:length(countgr), ranges = ranges,
            countgr = countgr, BPPARAM=BPPARAM,
            FUN = function(i, ranges, countgr){
                suppressPackageStartupMessages(require(FraseR))

                ## TODO init with NA since we dont extracted
                # we are missing counts for junctions spliced in other samples
                sample_count <- rep(0L, length(ranges))

                # get overlap and add counts to the corresponding ranges
                overlaps <- findOverlaps(countgr[[i]], ranges, type = "equal")
                sample_count[overlaps@to] <- mcols(countgr[[i]])$count

                return(sample_count)
            }
        )
    }

    # convert it to a DataFrame
    sample_counts_df <- DataFrame(
        matrix(unlist(sample_counts), ncol = length(sample_counts))
    )

    # merge it with the type columen and add it to the range object
    mcols(ranges) <- sample_counts_df

    # set correct naming
    colnames(mcols(ranges)) <- sample_names

    # return the object
    return(ranges)
}


##
## returns the name of the cache file if caching is enabled for the given sample
##
.getNonSplicedCountCacheFile <- function(sampleID, settings){
    # check if caching is enabled
    if(is.null(outputFolder(settings)) || outputFolder(settings) == "")
        return(NULL)

    # cache folder
    cachedir <- file.path(outputFolder(settings), "cache", "nonSplicedCounts")
    if(!dir.exists(cachedir)){
        dir.create(cachedir, recursive=TRUE)
    }

    # file name
    filename <- paste0("nonSplicedCounts-", sampleID, ".RDS")

    # return it
    return(file.path(cachedir, filename))
}

##
## counts non spliced reads based on the given target (acceptor/donor) regions
## TODO: 10k chunks hardcoded currently (needs some testing and a code to)
.countNonSplicedReads <- function(sampleID, spliceSiteCoords, settings,
                    internBPPARAM=SerialParam()){
    suppressPackageStartupMessages(library(FraseR))
    bamFile <- bamFiles(settings[samples(settings) == sampleID])[[1]]

    # check cache if available
    cacheFile <- .getNonSplicedCountCacheFile(sampleID, settings)
    if(!is.null(cacheFile) && file.exists(cacheFile)){
        # check if needs to be recalculated
        cache <- readRDS(cacheFile)
        numEqualOverlaps <- length(unique(to(findOverlaps(
                cache, spliceSiteCoords, type="equal"
        ))))
        if(numEqualOverlaps == length(spliceSiteCoords)){
            return(cache)
        } else {
            message(paste("The cache file does not contain the same",
                    "genomic positions. Remove cache and recounting ..."
            ))
            file.remove(cacheFile)
        }
    }

    # estimate chunk size
    rangeShift        <- 25*10^3
    numRangesPerChunk <- 100
    targetChunksTmp <- GenomicRanges::reduce(GenomicRanges::trim(suppressWarnings(
        GenomicRanges::shift(resize(spliceSiteCoords, width=rangeShift),
                    shift=-rangeShift/2
            )
    )))

    numChunks <- ceiling(max(1,length(targetChunksTmp)/numRangesPerChunk))
    chunkID <- rep(1:numChunks, each=numRangesPerChunk)[
                1:length(targetChunksTmp)
    ]
    targetChunks <- split(targetChunksTmp, chunkID)

    # extract the counts per chromosome
    countsList <- bplapply(targetChunks, bamFile=bamFile, settings=settings,
                spliceSites=spliceSiteCoords, BPPARAM=internBPPARAM,
                FUN=function(range, bamFile, settings, spliceSites){
# system.time({
        suppressPackageStartupMessages(library(FraseR))

        # restrict to the chromosome only
        param <- .mergeBamParams(scanBamParam(settings), range, TRUE)

        # extract raw data
        message(date(), ": Running on following chromosomes: '",
                paste(unique(seqnames(range)), collapse = ", "),
                "' \twith the number of ranges: '", length(range), "'",
                "\n\ton Bamfile: ",path(bamFile)
        )
        singleReadFrag <- readGAlignments(bamFile, param=param) %>%
                grglist() %>% reduce()

        # count data
        regionOfChunk <- subsetByOverlaps(spliceSites, range, type = "any")
        hits <- countOverlaps(regionOfChunk, singleReadFrag, minoverlap = 2)

        # clean memory
        #rm(regionOfChunk, singleReadFrag)
        #gc()
# }, TRUE)
        return(hits)
    })

    #
    mcols(spliceSiteCoords)$count <- unlist(countsList)
    spliceSiteCoords <- sort(spliceSiteCoords)

    # cache it if enabled
    if(!is.null(cacheFile)){
        saveRDS(spliceSiteCoords, cacheFile)
    }

    # return it
    return(spliceSiteCoords)
}


##
## extracts the splice site coordinates from a junctions GRange object
##
.extractSpliceSiteCoordinates <- function(junctions_gr, settings){
    if(settings@strandSpecific){
        splice_site_coords <- unlist(GRangesList(
            .extract_splice_site_coordinates_per_strand(junctions_gr, "+"),
            .extract_splice_site_coordinates_per_strand(junctions_gr, "-")
        ))
    } else {
        strand(junctions_gr) <- "*"
        splice_site_coords <-
            .extract_splice_site_coordinates_per_strand(junctions_gr, "*")
    }

    return(sort(unique(splice_site_coords)))
}


##
## extracts the splice site coordinates per strand from a
## given junctions GRange object
##
.extract_splice_site_coordinates_per_strand <- function(junctions_gr, strand){

    # get only the correct strand features
    junctions_gr <- junctions_gr[strand(junctions_gr) == strand,]

    # left side (acceptor on + and donor on -)
    left_side <- GRanges(
            seqnames = seqnames(junctions_gr),
            strand = strand(junctions_gr),
            ranges = IRanges(
                    start = start(junctions_gr) - 1,
                    end   = start(junctions_gr)
            ),
            seqlengths = seqlengths(junctions_gr),
            seqinfo = seqinfo(junctions_gr)
    )

    # right side (acceptor on - and donor on +)
    right_side <- GRanges(
            seqnames = seqnames(junctions_gr),
            strand = strand(junctions_gr),
            ranges = IRanges(
                    start = end(junctions_gr),
                    end   = end(junctions_gr) + 1
            ),
            seqlengths = seqlengths(junctions_gr),
            seqinfo = seqinfo(junctions_gr)
    )

    # annotate donor and acceptor sites
    if(strand == "+" | strand == "*"){
        mcols(left_side)$type = "Donor"
        mcols(right_side)$type = "Acceptor"
    } else {
        mcols(left_side)$type = "Acceptor"
        mcols(right_side)$type = "Donor"
    }

    return(sort(unlist(GRangesList(left_side, right_side))))
}

