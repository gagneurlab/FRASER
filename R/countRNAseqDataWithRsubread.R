#'
#' count with rsubread
#'

#'
#' creates a SAF data.table based on the given grange like object
#' @noRd
.GRange2SAF <- function(gr){
    data.table(
        GeneID  = seq_along(gr),
        Chr     = as.character(seqnames(gr)),
        Start   = start(gr),
        End     = end(gr),
        Strand  = as.character(strand(gr))
    )
}

#'
#' count non spliced reads with RSubRead
#'
.countNonSplicedReadsWithRsubread <- function(sampleID, settings, spliceSiteCoords, internBPPARAM){

    message(date(), ": Count non spliced reads for sample with Rsubread: ", sampleID)
    suppressPackageStartupMessages(library(FraseR))
    suppressPackageStartupMessages(library(Rsubread))

    bamFile <- bamFiles(settings[samples(settings) == sampleID])[[1]]

    # check cache if available
    cacheFile <- .getNonSplicedCountCacheFile(sampleID, settings)
    cacheFile <- gsub(".RDS$", "-rsubread.RDS", cacheFile)
    if(FALSE){#!is.null(cacheFile) && file.exists(cacheFile)){
        # check if needs to be recalculated
        cache <- try(readRDS(cacheFile), silent=TRUE)
        if(class(cache) == "try-error"){
            message(date(), ":Cache is currupt for sample:",
                    sampleID, ". Will recount it."
            )
            cache <- GRanges()
        }

        ov    <- findOverlaps(cache, spliceSiteCoords, type="equal")

        # we have all sites of interest cached
        if(all(1:length(spliceSiteCoords) %in% to(ov))){
            return(cache[from(ov)])
        } else {
            message(paste("The cache file does not contain the needed",
                      "genomic positions. Adding the remaining sites to",
                      "the cache ..."
        ))
        # get remaining splice sites
        spliceSiteCoords <-
            spliceSiteCoords[!(1:length(spliceSiteCoords) %in% to(ov))]
        }
    }

    internNcpu <- internBPPARAM
    if(is(internNcpu, "BiocParallelParam")){
        internNcpu <- bpworkers(internNcpu)
    }

    anno <- .GRange2SAF(spliceSiteCoords)
    res <- featureCounts(
        files = path(bamFile),
        annot.ext = anno,
        minOverlap = 2,

        allowMultiOverlap=TRUE,

        # multi-mapping reads
        countMultiMappingReads=TRUE,

        checkFragLength=FALSE,

        # read filtering
        minMQS=bamMapqFilter(scanBamParam(settings)),

        # strandness
        strandSpecific=settings@strandSpecific,

        # parameters specific to paired end reads
        isPairedEnd=TRUE,
        autosort=TRUE,
        nthreads=internNcpu
    )


    #
    mcols(spliceSiteCoords)$count <- res$counts[,1]
    spliceSiteCoords <- sort(spliceSiteCoords)


    # cache it if enabled (all splice sites!)
    if(!is.null(cacheFile)){
        if(exists("cache") && class(cache) == "GRanges"){
            message(date(), ": Adding splice sites to cache ...")
            saveRDS(unique(unlist(GRangesList(cache, spliceSiteCoords))),
                    cacheFile
            )
        } else {
            message(date(), ":Saving splice site cache ...")
            saveRDS(spliceSiteCoords, cacheFile)
        }
    }

    # merge with existing cache if needed
    if(exists("cache") && class(cache) == "GRanges"){
        # take only sites of interest
        cache <- cache[unique(from(ov))]

        # merge it with new sites
        spliceSiteCoords <- unlist(GRangesList(cache, spliceSiteCoords))
    }

    # return it
    return(spliceSiteCoords)
}
