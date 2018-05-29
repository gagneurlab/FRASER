########
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for annotating the ranges with
## biomaRt from ENSEMBL
##

#'
#' Annotates the given FraseRDataSet with the HGNC symbol with biomaRt
#'
#' @examples
#'
#' fds <- countRNAData(createTestFraseRSettings())
#' fds <- annotateRanges(fds)
#'
#' rowRanges(fds)[,"hgnc_symbol"]
#'
#' @export
annotateRanges <- function(fds, feature="hgnc_symbol", featureName=feature,
            biotype=list("protein_coding"), ensembl=NULL){

    # check input
    stopifnot(is(fds, "FraseRDataSet"))
    if(length(fds) == 0) return(fds)

    if(is.null(ensembl)){
        tryCatch({
            ensemblOutput <- capture.output(ensembl <- useEnsembl(
                    biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37
            ))
        },
        error=function(e){
                message("\nCheck if we have a internet connection!",
                        " Could not connect to ENSEMBL."
                )
        })
        if(is.null(ensembl)){
            message("Nothing was annotated!")
            return(fds)
        }
    }

    # check if USCS naming scheme is used
    useUSCS <- all(startsWith(seqlevels(fds), "chr"))

    # get gene annotation
    annotation <- getFeatureAsGRange(ensembl, feature, featureName,
            biotype, useUSCS)

    # annotate split reads
    mcols(fds, type="psi3")[[featureName]] <- getAnnotationFeature(
            data=rowRanges(fds), featureName, annotation
    )

    # annotate splice sites
    mcols(fds, type="psiSite")[[featureName]] <- getAnnotationFeature(
            data=rowRanges(nonSplicedReads(fds)),
            featureName, annotation
    )

    return(fds)
}


#'
#' use biomart to extract the current feature annotation
#'
#' @noRd
getFeatureAsGRange <- function(ensembl, feature, featureName,
                    biotype, useUSCS=TRUE){

    # contact biomaRt to retrive hgnc symbols
    ensemblResults <- getBM(
        attributes=c(feature, "chromosome_name",
                "start_position", "end_position"
        ),
        filters=c("biotype"),
        values=list(biotype=biotype),
        mart=ensembl
    )
    setnames(ensemblResults, feature, featureName)

    # remove emtpy symbols or non standard chromosomes
    ensemblResults <- ensemblResults[!is.na(ensemblResults[[featureName]]),]
    ensemblResults <- ensemblResults[ensemblResults[[featureName]] != "",]
    ensemblResults <- ensemblResults[
            !grepl("_|\\.", ensemblResults$chromosome_name),
    ]

    # create a grange object
    results <- makeGRangesFromDataFrame(
            ensemblResults, start.field="start_position",
            end.field="end_position", keep.extra.columns = TRUE
    )

    if(useUSCS){
        seqlevels(results) <- paste0("chr", seqlevels(results))
    }

    return(results)
}


#'
#' merge the retrieved annotations from a biomart with the ranges of our data
#'
#' @noRd
getAnnotationFeature <- function(data, feature, annotation){

    # find overlap with the query
    # suppress seqnames not in anothers object!
    # TODO: is there a nicer way to do this?
    suppressWarnings(hits <- findOverlaps(data, annotation))

    # extract only the feature and group them with a ";"
    featureDT <- data.table(
            from=from(hits),
            feature=mcols(annotation[to(hits)])[[feature]]
    )
    missingValues <- setdiff(1:length(data), unique(from(hits)))
    if(length(missingValues) > 0){
        featureDT <- rbind(
                featureDT, data.table(from=missingValues, feature=NA)
        )
    }

    # TODO only take first hit since it is tooo slow to merge them
    # with more then 10k entries
    featureDT <- unique(featureDT, by="from")
    # featureDT <- featureDT[,
    #        list(feature=paste(sort(unique(feature)), collapse = ";")),
    #        by="from"
    # ]

    return(featureDT[order(from),feature])
}


