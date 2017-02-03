########
## @author Christian Mertes \email{mertes@@in.tum.de}
## 
## This file contains all functions for annotating the ranges with 
## biomaRt from ENSEMBL
##

#'
#' Annotates the given FraseRDataSet with the HGNC symbol
#' with biomaRt
#' 
#' @export
annotateRanges <- function(dataset, feature="hgnc_symbol", biotype=list("protein_coding"), 
            ensembl=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)){
    
    # check if USCS naming scheme is used
    useUSCS <- all(startsWith(seqlevels(dataset@splitReads), "chr"))
    
    # get gene annotation
    annotation <- .getAllHgnsSymbolsAsGRange(ensembl, feature, biotype, useUSCS)
    
    # annotate split reads
    mcols(dataset@splitReads)[[feature]] <- .getAnnotationFeature(
            rowRanges(dataset@splitReads), feature, annotation
    )
    
    # annotate splice sites
    mcols(dataset@nonSplicedReads)[[feature]] <- .getAnnotationFeature(
            data=rowRanges(dataset@nonSplicedReads), 
            feature, annotation
    )
    
    return(dataset)    
}


#'
#' use biomart to extract the current HGNC annotation
#'
#' @noRd
.getAllHgnsSymbolsAsGRange <- function(ensembl, feature, biotype, useUSCS = TRUE){
    
    # contact biomaRt to retrive hgnc symbols
    ensemblResults <- getBM(
        attributes = c(feature, "chromosome_name", "start_position", "end_position"),
        filters = c("biotype"),
        values = list(biotype=biotype),
        mart = ensembl
    )
    
    # remove emtpy symbols or non standard chromosomes
    ensemblResults <- ensemblResults[!is.na(ensemblResults[[feature]]),]
    ensemblResults <- ensemblResults[ensemblResults[[feature]] != "",]
    ensemblResults <- ensemblResults[!grepl("_|\\.", ensemblResults$chromosome_name),]
    
    # create a grange object
    results <- makeGRangesFromDataFrame(ensemblResults, start.field = "start_position", 
                                        end.field = "end_position", keep.extra.columns = TRUE
    )
    
    if(useUSCS){
        seqlevels(results) <- paste0("chr", seqlevels(results))
    }
    
    return(results)
}


#'
#' retriev annotations from a biomart call
#'
#' @noRd
.getAnnotationFeature <- function(data, feature, annotation){
    
    # find overlap with the query
    # suppress seqnames not in anothers object!
    # TODO: is there a nicer way to do this?
    suppressWarnings(hits <- findOverlaps(data, annotation))
    
    # extract only the symbols and group them with a ";"
    featureDT <- data.table(
            from=from(hits), 
            feature=mcols(annotation[to(hits)])[[feature]]
    )
    missingValues <- setdiff(1:length(data), unique(from(hits)))
    if(length(missingValues) > 0){
        featureDT <- rbind(featureDT, data.table(from=missingValues, feature=NA))
    }
    # featureDT <- featureDT[,list(feature=paste(sort(unique(feature)), collapse = ";")),by="from"]
    # TODO only take first hit since it is tooo slow to merge them with more then 10k entries
    featureDT <- unique(featureDT, by="from")
    
    return(featureDT[order(from),feature])
}


