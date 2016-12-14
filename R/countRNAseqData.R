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
    message("Start counting the split reads ...")
    countList <- bplapply(settings@sampleData[,bamFile], 
                          FUN=.countSplitReads, 
                          settings=settings,
                          BPPARAM=settings@parallel,
                          internBPPARAM=internBPPARAM
    )
    names(countList) <- settings@sampleData[,sampleID]
    counts <- .mergeCounts(countList, settings@parallel)

    # count the retained reads
    message("Start counting the non spliced reads ...")
    countList <- bplapply(settings@sampleData[,bamFile], 
                          FUN=.countNonSplicedReads, 
                          settings=settings,
                          targets=granges(counts),
                          BPPARAM=settings@parallel,
                          internBPPARAM=internBPPARAM
    )
    names(countList) <- settings@sampleData[,sampleID]
    site_counts <- .mergeCounts(countList, settings@parallel)
    mcols(site_counts)$type=factor(countList[[1]]$type, 
                    levels = c("Acceptor", "Donor")
    )
    
    # create summarized objects
    splitCounts <- SummarizedExperiment(
        assays=list(rawCounts=mcols(counts)[settings@sampleData[,sampleID]]),
        rowRanges=granges(counts)
    )
    nonSplicedCounts <- SummarizedExperiment(
        assays=list(rawCounts=mcols(site_counts)[settings@sampleData[,sampleID]]),
        rowRanges=granges(site_counts)
    )
    
    # return it
    return(
        FraseRDataSet(
            settings=settings,
            splitReads=splitCounts,
            nonSplicedReads=nonSplicedCounts
        )
    )
}

##
## extracts the chromosomes within the given bamFile
##
.extractChromosomes <- function(bamFile){
    names(scanBamHeader(path(bamFile))[[path(bamFile)]]$target)
}

##
## count all split reads in a bam file
##
.countSplitReads <- function(bamFile, settings, internBPPARAM){
    suppressPackageStartupMessages(library(FraseR))
    
    # parallelize over chromosomes
    chromosomes <- FraseR:::.extractChromosomes(bamFile)
    
    # extract the counts per chromosome
    countsList <- bplapply(chromosomes, FUN=FraseR:::.countSplitReadsPerChromosome,
                       bamFile=bamFile, 
                       settings=settings,
                       BPPARAM=internBPPARAM
    )
    
    # sort and merge the results befor returning
    return(sort(unlist(GRangesList(countsList))))
}

##
## counting the split reads per chromosome
##
.countSplitReadsPerChromosome <- function(chromosome, bamFile, settings){
    # restrict to the chromosome only
    param <- ScanBamParam(which=GRanges(
            seqnames=chromosome, 
            ranges=IRanges(0, 536870912)
    ))

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
    mcols(junctionsCounts)$count <- countOverlaps(junctionsCounts,
                                                  junctions,
                                                  type = 'equal'
    )

    # sort it and return the GRange object
    return(sort(junctionsCounts))
}
    

## 
## merge multiple samples into one SummarizedExperiment object
## 
## TODO: how to init non found junctions/splice sites (0L or NA)? 
##
.mergeCounts <- function(counts, BPPARAM=SerialParam()){
  
  # prepare range object
  sample_names <- names(counts)
  counts <- GRangesList(counts)
  ranges <- sort(unique(unlist(counts)))
  mcols(ranges)$count <- NULL
  names(ranges) <- NULL
  
  # merge each sample counts into the combined range object
  sample_counts <- bplapply(1:length(counts), ranges = ranges, 
                            counts = counts, BPPARAM=BPPARAM,
          FUN = function(i, ranges, counts){
                suppressPackageStartupMessages(require(FraseR))
            
                # get sample name
                sample_name <- names(counts)[i]
                
                ## TODO init with NA since we dont extracted
                # we are missing counts for junctions spliced in other samples
                sample_count <- rep(0L, length(ranges))
                
                # get overlap and add counts to the corresponding ranges
                overlaps <- findOverlaps(counts[[i]], ranges, type = "equal")
                sample_count[overlaps@to] <- mcols(counts[[i]])$count
                #mcols(ranges[overlaps@to,])[[sample_name]] <- mcols(counts[[i]])$count
                
                return(sample_count)
          }
  )
  
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
## counts non spliced reads based on the given target (acceptor/donor) regions
## TODO: 10k chunks hardcoded currently (needs some testing and a code to)
.countNonSplicedReads <- function(bamFile, targets, settings, internBPPARAM=SerialParam()){
    suppressPackageStartupMessages(library(FraseR))
    
    # extract donor and acceptor sites
    splice_site_coordinates <- FraseR:::.extract_splice_site_coordinates(targets, settings)
     
    # estimate chunk size 
    numRangesPerChunk <- 50000
    numChunks <- ceiling(length(splice_site_coordinates)/numRangesPerChunk)
    chunkID <- rep(1:numChunks, each=numRangesPerChunk)[1:length(splice_site_coordinates)]
    targetChunks <- split(splice_site_coordinates, chunkID)   
    
    # extract the counts per chromosome
    countsList <- bplapply(targetChunks, bamFile=bamFile, 
                           settings=settings,
                           BPPARAM=internBPPARAM,
                           FUN=function(range, bamFile, settings){
                               single_read_fragments <- readGAlignments(bamFile, 
                                                param=ScanBamParam(which = range)) %>% 
                                                grglist() %>% reduce()
                               
                               countOverlaps(range, single_read_fragments, minoverlap = 2)
                           }
    )
	mcols(splice_site_coordinates)$count <- unlist(countsList)
	
	return(sort(splice_site_coordinates))
}


##
## extracts the splice site coordinates from a junctions GRange object
## 
.extract_splice_site_coordinates <- function(junctions_gr, settings){
	if(settings@strandSpecific){
		splice_site_coords <- unlist(GRangesList(
		    FraseR:::.extract_splice_site_coordinates_per_strand(junctions_gr, "+"),
		    FraseR:::.extract_splice_site_coordinates_per_strand(junctions_gr, "-")
		))
	} else { 
		strand(junctions_gr) <- "*"
		splice_site_coords <- FraseR:::.extract_splice_site_coordinates_per_strand(junctions_gr, "*")
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

