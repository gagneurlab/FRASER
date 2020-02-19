##
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for reading in data
## especially aligned RNA sequencing data
##


#' @title Count RNA-seq data
#'  
#' @description The FraseR package provides multiple functions to extract and 
#'               count both split and non-spliced reads from bam files. 
#'               See Detail and Functions for more information.
#'  
#' @details 
#' The functions described in this file extract and count both the 
#' split and the non-spliced reads from bam files. 
#' 
#' \code{\link{countRNAData}} is the main function that takes care of all 
#' counting steps and returns a FraseRDataSet containing the counts for all 
#' samples in the fds. 
#' 
#' \code{\link{getSplitReadCountsForAllSamples}} counts split reads for all 
#' samples and \code{\link{getNonSplitReadCountsForAllSamples}} counts non 
#' split reads overlapping splice sites for all samples. 
#' \code{\link{addCountsToFraseRDataSet}} adds these counts to an existing fds.
#' 
#' \code{\link{countSplitReads}} calculates the split read counts for a single
#' sample. \code{\link{countNonSplicedReads}} counts the non split reads 
#' overlapping with splice sites for a single sample.
#' 
#' \code{\link{mergeCounts}} merges the counts from different samples into a 
#' single count object, where the counts for junctions that are not present in 
#' a sample are set to zero. 
#' 
#' @param fds A FraseRDataSet object
#' @param NcpuPerSample A BiocParallel param object or a positive integer
#'             to configure the parallel backend
#'             of the internal loop per sample
#' @param junctionMap A object or file containing a map of
#'             all junctions of interest across all samples
#' @param minAnchor Minimum overlap around the Donor/Acceptor for
#'             non spliced reads. Default to 5
#' @param recount if TRUE the cache is ignored and the bam file is recounted.
#' @param genome the genome
#' @param BPPARAM the BiocParallel parameters for the parallelization
#' @param countDir The directory in which the tsv containing the position and 
#'                 counts of the junctions should be placed.
#' @param countFiles If specified, the split read counts for all samples are 
#'                  read from the specified files. Should be a vector of paths 
#'                  to files containing the split read counts for the 
#'                  individual samples. Reading from files is only supported 
#'                  for tsv(.gz) or RDS files containing GRranges objects. The 
#'                  order of the individual sample files should correspond to 
#'                  the order of the samples in the fds.
#' @param outFile The full path to the output tsv file containing the merged 
#'               counts. If the given file already exists, this counts from 
#'               this file will be read in and used in the following (i.e. the 
#'               reads are not recounted), unless the option recount=TRUE is 
#'               used. If this file doesn't exist or if recount=TRUE, then it 
#'               will be created after counting has finished.
#' @param splitCountRanges The merged GRanges object containing the positions 
#'              of all the introns in the dataset over all samples.
#' @param splitCounts The SummarizedExperiment object containing the 
#'              position and counts of all the introns in the dataset 
#'              for all samples.
#' @param nonSplitCounts The SummarizedExperiment object containing the 
#'              position and non split read counts of all splice sites 
#'              present in the dataset for all samples.
#' @param sampleID The ID of the sample to be counted.
#' @param spliceSiteCoords A GRanges object containing the positions of the 
#'         splice sites. If it is NULL, then splice sites coordinates are 
#'         calculated first based on the positions of the junctions defined 
#'         from the split reads.
#' @param longRead If TRUE, then the isLongRead option of 
#'         Rsubread::featureCounts is used when counting the non spliced reads 
#'         overlapping splice sites.
#' @param countList A list of GRanges objects containing the counts that should
#'         be merged into one object.
#' @param assumeEqual Logical indicating whether all objects in 
#'         \code{countList} can be assumed to contain counts for the same ranges. 
#'         If FALSE, merging of the ranges is performed.
#' @param ... Further parameters passed on to Rsubread::featureCounts.
#' 
#' @name countRNA
#' @rdname countRNA
#' 
#' @examples
#'   fds <- countRNAData(createTestFraseRSettings())
#'
NULL


#' @describeIn countRNA This method extracts and counts the split reads and
#'             non spliced reads from RNA bam files.
#' @return \code{\link{countRNAData}} returns a FraseRDataSet.
#' 
#' @export
countRNAData <- function(fds, NcpuPerSample=1, junctionMap=NULL, minAnchor=5,
                        recount=FALSE, BPPARAM=bpparam(), genome=NULL,
                        countDir=file.path(workingDir(fds), "savedObjects", 
                                            nameNoSpace(name(fds))),
                        ...){
    
    # Check input TODO
    stopifnot(is(fds, "FraseRDataSet"))
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
    
    # count splitreads first for every sample
    splitCounts <- getSplitReadCountsForAllSamples(fds=fds, 
                                                NcpuPerSample=NcpuPerSample, 
                                                junctionMap=junctionMap, 
                                                recount=recount, 
                                                BPPARAM=BPPARAM, 
                                                genome=genome,
                                                outFile=file.path(countDir, 
                                                        "splitCounts.tsv.gz"))
    
    # count non spliced reads for every samples
    nonSplicedCounts <- getNonSplitReadCountsForAllSamples(fds=fds, 
                                                splitCountRanges=
                                                    rowRanges(splitCounts), 
                                                NcpuPerSample=NcpuPerSample, 
                                                minAnchor=minAnchor,
                                                recount=recount, 
                                                BPPARAM=BPPARAM,
                                                outFile=file.path(countDir, 
                                                    "nonSplitCounts.tsv.gz"),
                                                ...)
    
    # create final FraseR dataset
    fds <- addCountsToFraseRDataSet(fds, splitCounts, nonSplicedCounts)
    gc()
    
    # return it
    return(fds)
}

#' @describeIn countRNA This method creates a GRanges 
#'             object containing the split read counts from all 
#'             specified samples.  
#' @return \code{\link{getSplitReadCountsForAllSamples}} returns a GRanges 
#' object.
#' 
#' @export
getSplitReadCountsForAllSamples <- function(fds, NcpuPerSample=3, 
                                            junctionMap=NULL, recount=FALSE, 
                                            BPPARAM=bpparam(), genome=NULL, 
                                            countFiles=NULL,
                                            outFile=file.path(workingDir(fds), 
                                                        "savedObjects", 
                                                        nameNoSpace(name(fds)),
                                                        "splitCounts.tsv.gz")){
    
    # check if outFile with mergedCounts already exists
    # if so, don't recalculate the split counts
    if(file.exists(outFile) && isFALSE(recount)){
        
        counts <- makeGRangesFromDataFrame(fread(outFile), 
                                            keep.extra.columns=TRUE)
        
        # if counts for all samples are present in the tsv, return those counts
        if(all(samples(fds) %in% colnames(mcols(counts)))){
            message(date(), ": File with the merged split read counts exists ", 
                    "already and will be used. If you want to count the split ",
                    "reads again, use the option recount=TRUE or remove this ",
                    "file: \n", outFile)
            # create summarized object for counts
            h5 <- saveAsHDF5(fds, "rawCountsJ", 
                                as.matrix(mcols(counts)[samples(fds)]))
            colnames(h5) <- samples(fds)
            final_splitCounts <- SummarizedExperiment(
                colData=colData(fds),
                assays=list(rawCountsJ=h5),
                rowRanges=counts[,!colnames(mcols(counts)) %in% samples(fds)]
            )
            rm(h5)
            gc()
            return(final_splitCounts)
        }
    }
    
    # split reads need to be counted for the given samples
    if(is.null(countFiles)){ 
        
        # count splitreads first for every sample
        message(date(), ": Start counting the split reads ...")
        countList <- bplapply(samples(fds),
                                FUN=countSplitReads,
                                fds=fds,
                                BPPARAM=BPPARAM,
                                NcpuPerSample=NcpuPerSample,
                                genome=genome,
                                recount=recount)
        
        
    } else { # counts should be read from files
        
        # check that all given files exist and that there is a file for every
        # sample
        stopifnot(all(vapply(countFiles, file.exists, FUN.VALUE = logical(1))))
        stopifnot(length(samples(fds)) == length(countFiles))
        
        message(date(), ": Reading in the files with the split reads ...")
            
        if(all(grepl(x=countFiles, pattern=".tsv"))){
                
            # read in tsv files as GRanges objects
            countList <- lapply(countFiles, 
                                function(x){
                                    makeGRangesFromDataFrame(
                                        fread(x), keep.extra.columns=TRUE) 
                                } )
            
        } else if(all(grepl(x=countFiles, pattern=".RDS", ignore.case=TRUE))){
                
            # read in rds files with GRanges objects
            countList <- lapply(countFiles, readRDS)
            stopifnot(all(vapply(countList, FUN=is, class2="GRanges", 
                                FUN.VALUE = logical(1))))
                
        } else{
            stop(paste0("Reading from files is only supported for tsv(.gz) or ",
                        "RDS files containing GRranges objects."))
        }
        
    }
    
    names(countList) <- samples(fds)
    BPPARAM_old <- setMaxThreads(BPPARAM, 10, 10)
    counts <- mergeCounts(countList, fds, junctionMap=junctionMap,
                            assumeEqual=FALSE, BPPARAM=BPPARAM_old$BPPARAM)
    BPPARAM <- setMaxThreads(BPPARAM_old$BPPARAM,
                            BPPARAM_old$numWorkers, BPPARAM_old$numTasks, 
                            set=TRUE)
    
    rm(countList)
    gc()
    
    # splice site map
    rowRanges(counts) <- annotateSpliceSite(rowRanges(counts))
    
    # save summarized experiment
    outDir <- file.path(dirname(outFile), "splitCounts")
    if(!dir.exists(outDir)){
        dir.create(outDir, recursive=TRUE)
    }
    message(date(), ": Writing split counts to folder: ", outDir)
    saveHDF5SummarizedExperiment(counts, dir=outDir, replace=TRUE)
    
    # write tsv
    # writeCountsToTsv(assay(counts), file=outFile)
    
    return(counts)
}

#' @describeIn countRNA This method creates a GRanges 
#'              object containing the non split read counts at the 
#'              exon-intron boundaries inferred from the GRanges object 
#'              containing the positions of all the introns in this dataset.  
#' @return \code{\link{getNonSplitReadCountsForAllSamples}} returns a 
#'          GRanges object.
#' @export
getNonSplitReadCountsForAllSamples <- function(fds, splitCountRanges, 
                    NcpuPerSample=3, minAnchor=5, recount=FALSE, 
                    BPPARAM=bpparam(), longRead=TRUE, outFile=file.path(
                            workingDir(fds), "savedObjects", 
                            nameNoSpace(name(fds)), "nonSplitCounts.tsv.gz")){
    
    # check if outFile with mergedCounts already exists
    # if so, don't recalculate the non split counts
    if(file.exists(outFile) && isFALSE(recount)){
        siteCounts <- makeGRangesFromDataFrame(fread(outFile), 
                keep.extra.columns=TRUE)
        
        # if counts for all samples are present in the tsv, return those counts
        if(all(samples(fds) %in% colnames(mcols(siteCounts)))){
            message(date(), ": File with the merged non-split read counts ", 
                    "exists already and will be used. If you want to count ",
                    "the non-split reads again, use the option recount=TRUE ",
                    "or remove this file: \n", outFile)
            h5 <- saveAsHDF5(fds, "rawCountsSS",
                                as.matrix(mcols(siteCounts)[samples(fds)]))
            colnames(h5) <- samples(fds)
            final_nonSplicedCounts <- SummarizedExperiment(
                colData=colData(fds),
                assays=list(rawCountsSS=h5),
                rowRanges=siteCounts[,!colnames(mcols(siteCounts)) %in% 
                                            samples(fds)]
            )
            rm(h5)
            gc()
            return(final_nonSplicedCounts)
        }
    }
    
    if(!("startID" %in% colnames(mcols(splitCountRanges))) | 
        !("endID" %in% colnames(mcols(splitCountRanges)))){
        # splice site map
        splitCountRanges <- annotateSpliceSite(splitCountRanges)
    }
    
    
    # count the retained reads
    message(date(), ": Start counting the non spliced reads ...")
    message(date(), ": In total ", length(splitCountRanges),
            " splice junctions are found.")
    
    # extract donor and acceptor sites
    spliceSiteCoords <- extractSpliceSiteCoordinates(splitCountRanges, fds)
    message(date(), ": In total ", length(spliceSiteCoords),
            " splice sites (acceptor/donor) will be counted ...")
    
    # count non spliced reads for every samples
    saveSpliceSiteCoordinates(spliceSiteCoords, fds)
    countList <- bplapply(samples(fds),
                            FUN=countNonSplicedReads,
                            fds=fds,
                            splitCountRanges=splitCountRanges,
                            BPPARAM=BPPARAM,
                            NcpuPerSample=NcpuPerSample,
                            minAnchor=minAnchor,
                            recount=recount,
                            spliceSiteCoords=spliceSiteCoords,
                            longRead=longRead)
    names(countList) <- samples(fds)
    siteCounts <- mergeCounts(countList, fds=fds, assumeEqual=TRUE, 
                                spliceSiteCoords=spliceSiteCoords )
    
    rm(countList)
    gc()
    
    # save summarized experiment
    outDir <- file.path(dirname(outFile), "nonSplitCounts")
    if(!dir.exists(outDir)){
        dir.create(outDir, recursive=TRUE)
    }
    message(date(), ": Writing non-split counts to folder: ", outDir)
    saveHDF5SummarizedExperiment(siteCounts, dir=outDir, replace=TRUE)
    
    # write tsv
    # writeCountsToTsv(assay(siteCounts), file=outFile)
    
    return(siteCounts)
    
}


#' @describeIn countRNA This method adds the split read and 
#'              non split read counts to a existing FraseRDataSet 
#'              containing the settings.  
#' @return \code{\link{addCountsToFraseRDataSet}} returns a FraseRDataSet.
#' @export
addCountsToFraseRDataSet <- function(fds, splitCounts, nonSplitCounts){
    
    # create final FraseR dataset
    fds <- new("FraseRDataSet",
                splitCounts,
                name            = name(fds),
                bamParam        = scanBamParam(fds),
                strandSpecific  = strandSpecific(fds),
                workingDir      = workingDir(fds),
                nonSplicedReads = nonSplitCounts,
                metadata        = metadata(fds)
    )
    
    # save it so the FraseR object also gets saved
    fds <- saveFraseRDataSet(fds)
    
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


#' @describeIn countRNA This method counts all split reads in a 
#'     bam file for a single sample.
#' @return \code{\link{countSplitReads}} returns a GRanges object.
#' @export
countSplitReads <- function(sampleID, fds, NcpuPerSample=1, genome=NULL, 
                            recount=FALSE){
    bamFile <- bamFile(fds[,samples(fds) == sampleID])[[1]]
    
    # check cache if available
    cacheFile <- getSplitCountCacheFile(sampleID, fds)
    if(isFALSE(recount) && !is.null(cacheFile) && file.exists(cacheFile)){
        cache <- readRDS(cacheFile)
        bamWhich <- unlist(bamWhich(scanBamParam(fds)))
        if(length(bamWhich) > 0){
            userTargetGR <- GRanges(seqnames=names(unlist(bamWhich)),
                                    ranges=unlist(bamWhich), strand="*")
            from <- unique(from(findOverlaps(cache, userTargetGR)))
            cache <- cache[from]
        }
        if(length(cache) > 0){
            message(date(), ": Using existing split read counts for sample: ", 
                    sampleID)
            return(checkSeqLevelStyle(cache, fds, sampleID,
                                        sampleSpecific=FALSE))
        }
    }
    message(date(), ": Count split reads for sample: ", sampleID)
    
    # parallelize over chromosomes
    chromosomes <- extractChromosomes(bamFile)
    
    if(length(genome) > 1){
        genome <- genome[sampleID]
    }
    
    # extract the counts per chromosome
    countsList <- bplapply(chromosomes, FUN=countSplitReadsPerChromosome,
            bamFile=bamFile, settings=fds, genome=genome,
            BPPARAM=getBPParam(NcpuPerSample, length(chromosomes)))
    
    # sort and merge the results befor returning/saving
    countsGR <- sort(unlist(GRangesList(countsList)))
    saveRDS(countsGR, cacheFile)
    
    return(checkSeqLevelStyle(countsGR, fds, sampleID, sampleSpecific=FALSE))
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


#' @describeIn countRNA This method merges counts for multiple 
#'                      samples into one SummarizedExperiment object.
#' @return \code{\link{mergeCounts}} returns a SummarizedExperiment object.
#' @export
mergeCounts <- function(countList, fds, junctionMap=NULL, assumeEqual=FALSE,
                        spliceSiteCoords=NULL, BPPARAM=SerialParam()){
    
    # prepare range object
    sample_names <- names(countList)
    
    if(assumeEqual){
        stopifnot(is(spliceSiteCoords, "GRanges"))
        ranges <- spliceSiteCoords
        mcols(ranges)$count <- NULL
        names(ranges) <- NULL
        mcols(ranges)$type <- factor(ranges$type, levels=c("Acceptor", "Donor"))
        
        message(paste(date(), ": Fast merging of counts ..."))
        sample_counts <- countList
    } else {
        
        if(!"SeqLevelStyle" %in% colnames(colData(fds))){
            colData(fds)[,"SeqLevelStyle"] <- 
                vapply(bamFile(fds), FUN.VALUE=character(1), 
                       FUN=function(bamFile){
                           seqlevelsStyle(BamFile(bamFile, yieldSize = 2e6))
                        }
                       )
        }
        countList <- mapply(countList, samples(fds), 
                            FUN=function(gr, sampleID){
                                checkSeqLevelStyle(gr, fds, sampleID, 
                                                   sampleSpecific=FALSE)
                                }
                            )
        
        countList <- uniformSeqInfo(countList)
        countgr <- GRangesList(countList)
        if(!is.null(junctionMap)){
            stopifnot(is(junctionMap, "GRanges"))
            countgr <- c(countgr, GRangesList(junctionMap))
        }
        ranges <- sort(unique(unlist(countgr)))
        mcols(ranges)$count <- NULL
        names(ranges) <- NULL
        
        message(paste(date(), ": count ranges need to be merged ..."))
        # merge each sample counts into the combined range object
        sample_counts <- bplapply(countList, ranges = ranges, BPPARAM=BPPARAM,
                FUN = function(gr, ranges){
                        
                        # init with 0 since we did not find any read
                        # for this sample supporting a given site
                        sample_count <- integer(length(ranges))
                        
                        # get overlap and add counts to the 
                        # corresponding ranges
                        overlaps <- findOverlaps(gr, ranges, type="equal")
                        sample_count[overlaps@to] <- mcols(gr)$count
                        
                        return(sample_count)
                }
        )
    }
    
    # merge it with the type column and add it to the range object
    counts <- do.call(cbind, sample_counts)
    colnames(counts) <- sample_names
    
    # create assay for summarized object
    aName <- ifelse(assumeEqual, "rawCountsSS", "rawCountsJ")
    h5 <- saveAsHDF5(fds, aName, as.matrix(counts[,samples(fds)]))
    colnames(h5) <- samples(fds)
    aList <- list(a=h5)
    names(aList) <- aName
    
    # return the object
    final_counts <- SummarizedExperiment(
            colData=colData(fds), assays=aList, rowRanges=ranges)
    return(final_counts)
}


#'
#' returns the name of the cache file if caching is enabled for the given sample
#' @noRd
getNonSplicedCountCacheFile <- function(sampleID, settings){
    
    # cache folder
    cachedir <- getNonSplicedCountCacheFolder(settings)
    
    # file name
    filename <- paste0("nonSplicedCounts-", sampleID, ".h5")
    
    # return it
    return(file.path(cachedir, filename))
}

#' 
#' Save splice site coordinates (GRanges) as RDS
#' @noRd
saveSpliceSiteCoordinates <- function(spliceSiteCoords, settings){
    # cache folder
    cachedir <- getNonSplicedCountCacheFolder(settings)
    # file name
    filename <- "spliceSiteCoordinates.RDS"
    
    # save it
    saveRDS(sort(spliceSiteCoords), file.path(cachedir, filename))
}

#'
#' returns the name of the cache folder if caching is enabled 
#' @noRd
getNonSplicedCountCacheFolder <- function(settings){
    
    # cache folder
    cachedir <- file.path(workingDir(settings), "cache", "nonSplicedCounts", 
                            nameNoSpace(name(settings)))
    if(!dir.exists(cachedir)){
        dir.create(cachedir, recursive=TRUE)
    }
    
    # return it
    return(cachedir)
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



#' @describeIn countRNA This method counts non spliced reads based 
#'     on the given target (acceptor/donor) regions for a single sample.
#' @return \code{\link{countNonSplicedReads}} returns a GRanges object.
#' @export
countNonSplicedReads <- function(sampleID, splitCountRanges, fds,
                    NcpuPerSample=3, minAnchor=5, recount=FALSE,
                    spliceSiteCoords=NULL, longRead=TRUE){
    
    if(is.null(spliceSiteCoords) | !is(spliceSiteCoords, "GRanges")){
        
        # splice site map
        if(!("startID" %in% colnames(mcols(splitCountRanges))) | 
            !("endID" %in% colnames(mcols(splitCountRanges)))){
            splitCountRanges <- annotateSpliceSite(splitCountRanges)
        }
        
        # extract donor and acceptor sites
        spliceSiteCoords <- extractSpliceSiteCoordinates(splitCountRanges, fds)
    }
    
    
    bamFile <- bamFile(fds[,samples(fds) == sampleID])[[1]]
    
    # check cache if available
    cacheFile <- getNonSplicedCountCacheFile(sampleID, fds)
    if(isFALSE(recount) && !is.null(cacheFile) && file.exists(cacheFile)){
        # check if needs to be recalculated
        cache <- try(HDF5Array(cacheFile, "nonSplicedCounts"), silent=TRUE)
        if(is(cache, "try-error")){
            message(date(), ":Cache is currupt for sample:",
                    sampleID, ". Will recount it."
            )
            cache <- matrix()
        }
        # cache <- checkSeqLevelStyle(cache, fds, sampleID, FALSE)
        # ov    <- findOverlaps(cache, spliceSiteCoords, type="equal")
        
        # we have all sites of interest cached
        if(length(spliceSiteCoords) == nrow(cache)){
            message(date(), 
                    ": Using existing non spliced read counts for sample: ", 
                    sampleID)
            return(cache)
        } else {
            message(paste("The cache file does not contain the needed",
                            "genomic positions. Adding the remaining sites to",
                            "the cache ..."
            ))
        }
    }
    
    message(date(), ": Count non spliced reads for sample: ", sampleID)
    
    # extract the counts with Rsubread
    tmp_ssc <- checkSeqLevelStyle(spliceSiteCoords, fds, sampleID, TRUE)
    anno <- GRanges2SAF(tmp_ssc, minAnchor=minAnchor)
    rsubreadCounts <- featureCounts(files=bamFile, annot.ext=anno,
            minOverlap=minAnchor*2, 
            allowMultiOverlap=TRUE,
            checkFragLength=FALSE,
            minMQS=bamMapqFilter(scanBamParam(fds)),
            strandSpecific=as.integer(strandSpecific(fds)),
            
            # activating long read mode
            isLongRead=longRead,
            
            # multi-mapping reads
            countMultiMappingReads=TRUE,
            
            # for counting only non spliced reads we 
            # skip this information
            isPairedEnd=FALSE,
            
            # this is important, otherwise it will sort it by name
            autosort=FALSE,
            nthreads=NcpuPerSample,
            tmpDir=file.path(file_path_as_absolute(workingDir(fds)), "cache")
    )
    
    # extract results
    mcols(spliceSiteCoords)$count <- rsubreadCounts$counts[,1]
    spliceSiteCoords <- sort(spliceSiteCoords)
    
    # get counts that will be cached
    cache <- as.matrix(ncol=1, spliceSiteCoords$count)
    rowChunkSize <- min(nrow(cache), options()[['FraseR-hdf5-chunk-nrow']])
    
    # cache counts as hdf5 file
    message("Saving splice site cache ...")
    if(file.exists(cacheFile)){
        unlink(cacheFile)
    }
    writeHDF5Array(cache, filepath=cacheFile, name="nonSplicedCounts", 
                    chunkdim=c(rowChunkSize,1), level=7, verbose=FALSE)
    
    # get counts as DelayedMatrix
    sample_counts <- HDF5Array(filepath=cacheFile, name="nonSplicedCounts")
        
    # return it
    return(sample_counts)
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
            stopifnot(is(map, "GRanges"))
            return(map)
        }
        message(date(), ": currently not supported.")
        message(date(), ": will only use the provided junctions")
        return(NULL)
    }
    if(is(junctionMap, "GRanges")){
        return(junctionMap)
    }
    stop("Objecttype (class) for junction map currently not supported:",
        class(junctionMap)
    )
}

#' extracts the splice site coordinates from a junctions GRange object (
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

#' annotates the given GRange object with unique identifier for splice sites 
#' (acceptors and donors).
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
    annotadedDT[,id:=seq_len(nrow(annotadedDT))]
    
    # convert back to granges
    annogr <- makeGRangesFromDataFrame(annotadedDT, keep.extra.columns=TRUE)
    
    ids <- lapply(c("start", "end"), function(type){
        # reduce annogr to only the specific type to prevent overlap
        annogrtmp <- annogr[annogr$type == type]
        ov <- findOverlaps(gr, annogrtmp, type=type)
        mcols(annogrtmp[to(ov)])[["id"]]
    })
    
    names(ids) <- paste0(c("start", "end"), "ID")
    mcols(gr)[names(ids)] <- DataFrame(ids)
    
    return(gr)
}

setMaxThreads <- function(BPPARAM, maxWorkers=bpworkers(BPPARAM), 
                            maxTasks=bptasks(BPPARAM), set=FALSE){
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

#' writes a GRanges object with the counts as a tsv (or tsv.gz) file.
#' @noRd
writeCountsToTsv <- function(counts, file="counts.tsv.gz"){
    if(!dir.exists(dirname(file))){
        dir.create(dirname(file), recursive=TRUE)
    }
    message(date(), ": Writing counts to file: ", file)
    fwrite(as.data.table(counts), file=file, sep = "\t")
}
