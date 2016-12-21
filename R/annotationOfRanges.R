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
annotateRanges <- function(dataset, biotype = list("protein_coding"), 
            ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)){
    
    # annotate split reads
    mcols(dataset@splitReads)$hgnc_symbol      <- .getHgncSymbols(dataset@splitReads, biotype, ensembl)
    
    # annotate splice sites
    mcols(dataset@nonSplicedReads)$hgnc_symbol <- .getHgncSymbols(dataset@nonSplicedReads, biotype, ensembl)

    return(dataset)    
}

#'
#'
#' @noRd
.getHgncSymbols <- function(ranges, biotype, ensembl){
    
    # get chormosomal regions
    regionsAsString <- paste0(seqnames(ranges), ":", start(ranges), "-", end(ranges))
    if(all(startsWith(seqlevels(ranges), "chr"))){
        regionsAsString <- gsub("^chr", "", regionsAsString)
    }
    
    # contact biomaRt to retrive hgnc symbols
    ensemblResults <- getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
        filters = c("chromosomal_region", "biotype"),
        values = list(chromosomal_region=regionsAsString, biotype=biotype),
        mart = ensembl
    )
    
    # remove emtpy symbols or non standard chromosomes
    ensemblResults <- ensemblResults[ensemblResults$hgnc_symbol != "",]
    ensemblResults <- ensemblResults[!grepl("_|\\.", ensemblResults$chromosome_name),]
    
    # create a grange object
    results <- makeGRangesFromDataFrame(ensemblResults, start.field = "start_position", 
            end.field = "end_position", keep.extra.columns = TRUE
    )
    if(all(startsWith(seqlevels(ranges), "chr"))){
        seqlevels(results) <- paste0("chr", seqlevels(results))
    }
    
    # find overlap with the query
    hits <- findOverlaps(ranges, results) 
    
    # extract only the symbols and group them with a ";"
    hgnc_symbols <- data.table(from=from(hits), hgnc_symbol=results[to(hits)]$hgnc_symbol)
    missing_values <- setdiff(1:length(ranges), unique(from(hits)))
    if(length(missing_values) > 0){
        hgnc_symbols <- rbind(hgnc_symbols, data.table(from=missing_values, hgnc_symbol=NA))
    }
    hgnc_symbols <- hgnc_symbols[,list(symbols=paste(sort(unique(hgnc_symbol)), collapse = ";")),by="from"][,symbols]
    
    return(hgnc_symbols)
}

#'
#' use biomart to extract the current HGNC annotation
#'
#' @noRd
.getAllHgnsSymbolsAsGRange <- function(ensembl, biotype, useUCSCScheme = TRUE){
    
    # contact biomaRt to retrive hgnc symbols
    ensemblResults <- getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
        filters = c("biotype"),
        values = list(biotype=biotype),
        mart = ensembl
    )
    
    # remove emtpy symbols or non standard chromosomes
    ensemblResults <- ensemblResults[ensemblResults$hgnc_symbol != "",]
    ensemblResults <- ensemblResults[!grepl("_|\\.", ensemblResults$chromosome_name),]
    
    # create a grange object
    results <- makeGRangesFromDataFrame(ensemblResults, start.field = "start_position", 
                                        end.field = "end_position", keep.extra.columns = TRUE
    )
    
    if(useUCSCScheme){
        seqlevels(results) <- paste0("chr", seqlevels(results))
    }
    
    return(results)
}


