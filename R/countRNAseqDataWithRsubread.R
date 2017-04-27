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

.countNonSplicedReadsWithRsubread <- function(sampleID, settings, spliceSiteCoords, internBPPARAM){

    .GRange2SAF <- function(gr){
        data.table(
            GeneID  = seq_along(gr),
            Chr     = as.character(seqnames(gr)),
            Start   = start(gr),
            End     = end(gr),
            Strand  = as.character(strand(gr))
        )
    }

    message(date(), ": Count non spliced reads for sample with Rsubread: ", sampleID)
    suppressPackageStartupMessages(library(FraseR))
    suppressPackageStartupMessages(library(Rsubread))

    bamFile <- bamFiles(settings[samples(settings) == sampleID])[[1]]

    # check cache if available
    cacheFile <- .getNonSplicedCountCacheFile(sampleID, settings)
    cacheFile <- gsub(".RDS$", "-rsubread.RDS", cacheFile)
    if(!is.null(cacheFile) && file.exists(cacheFile)){
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
    tmp_dir <- file.path(outputFolder(settings), "cache", "rsubread-tmp")
    if(!dir.exists(tmp_dir)) {
        dir.create(tmp_dir, recursive=TRUE)
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
        strandSpecific=0,

        # parameters specific to paired end reads
        isPairedEnd=FALSE,
        autosort=TRUE,
        nthreads=internNcpu,
        tmpDir=tmp_dir
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
