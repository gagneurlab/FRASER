########
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for annotating the ranges with
## biomaRt from ENSEMBL
##

#'
#' Annotates the given FraseRDataSet with the HGNC symbol with biomaRt
#' 
#' @return FraseRDataSet
#' 
#' @examples
#'
#' fds <- countRNAData(createTestFraseRSettings())
#' fds <- annotateRanges(fds)
#'
#' rowRanges(fds, type="psi5")[,"hgnc_symbol"]
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
    for(i in c("psi3", "psiSite")){
        gr <- rowRanges(fds, type=i)
        if(any(strand(gr) == "*")){
            strand(annotation) <- "*"
        }
        annos <- getAnnotationFeature(data=gr, featureName, annotation)
        mcols(fds, type=i)[[featureName]] <- annos
    }

    return(fds)
}


#'
#' use biomart to extract the current feature annotation
#'
#' @noRd
getFeatureAsGRange <- function(ensembl, feature, featureName,
                    biotype, useUSCS=TRUE){

    # contact biomaRt to retrive hgnc symbols
    ensemblResults <- getBM(attributes=c(feature, "chromosome_name",
                    "start_position", "end_position", "strand"),
            filters=c("biotype"), values=c(biotype), mart=ensembl)
    setnames(ensemblResults, feature, featureName)

    # remove emtpy symbols or non standard chromosomes
    ensemblResults <- ensemblResults[!is.na(ensemblResults[[featureName]]),]
    ensemblResults <- ensemblResults[ensemblResults[[featureName]] != "",]
    ensemblResults <- ensemblResults[
            !grepl("_|\\.", ensemblResults$chromosome_name),]
    ensemblResults[,"strand"]  <- ifelse(
            ensemblResults[,"strand"] < 0, "-", "+")
    
    # create a grange object
    results <- makeGRangesFromDataFrame(ensemblResults,
            start.field="start_position", end.field="end_position", 
            strand.field="strand", keep.extra.columns=TRUE)

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
    missingValues <- setdiff(seq_along(data), unique(from(hits)))
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


#'
#' Find annotated junctions
#' 
#' Annotate the object with a given annotation. 
#' 
#' @return FraseRDataSet
#' 
#' @export
findAnnotatedJunction <- function(fds, annotation, annotateNames=TRUE,
                    seqLevelStyle=seqlevelsStyle(fds)[1],
                    stranded=strandSpecific(fds), ...){
    
    # import annotation if given as file
    if(isScalarCharacter(annotation)){
        txdb <- makeTxDbFromGFF(annotation, ...)
    }
    
    # extract introns
    intronsBT <- intronsByTranscript(txdb, use.names=TRUE)
    introns <- unlist(intronsBT)
    
    # check if strandspecific data is used
    gr <- rowRanges(fds, type="psi5")
    
    if(isFALSE(stranded)){
        strand(gr) <- "*"
    }
    
    # set correct seq level names
    seqlevelsStyle(introns) <- seqLevelStyle
    seqlevelsStyle(gr) <- seqLevelStyle
    
    # find overlaps with introns
    ov <- findOverlaps(gr, introns, type="equal")
    newAnno <- data.table(idx=seq_row(fds), known_intron=FALSE, 
            known_start=FALSE, known_end=FALSE, TranscriptNames=NA_character_)
    newAnno[idx %in% from(ov), known_intron:=TRUE]
    
    if(isTRUE(annotateNames)){
        ovanno <- as.data.table(ov)
        ovanno[,name:=names(introns)[subjectHits]]
        ovanno <- ovanno[,.(
                TranscriptNames=paste(unique(name), collapse=",")),
                by=queryHits]
        newAnno <- merge(newAnno, ovanno, by.x="idx", 
                by.y="queryHits", all.x=TRUE)
    }
    
    # overlap start/stop
    ov <- findOverlaps(gr, introns, type="start")
    ovdt <- data.table(from=from(ov), to=to(ov), 
            strand=as.vector(strand(introns)[to(ov)]))[,
                    .(strand=strand[1]),by=from]
    newAnno[ovdt$from, known_start:=known_start | ovdt$strand %in% c("+", "*")]
    newAnno[ovdt$from, known_end  :=known_end   | ovdt$strand %in% c("-")]
    
    ov <- findOverlaps(gr, introns, type="end")
    ovdt <- data.table(from=from(ov), to=to(ov), 
            strand=as.vector(strand(introns)[to(ov)]))[,
                    .(strand=strand[1]),by=from]
    newAnno[ovdt$from, known_start:=known_start | ovdt$strand %in% c("-")]
    newAnno[ovdt$from, known_end  :=known_end   | ovdt$strand %in% c("+", "*")]
    
    # update fds object
    columns2Write <- c("known_intron", "known_start", "known_end", 
            "TranscriptNames")
    col2Take <- !colnames(mcols(fds, type="psi5")) %in% columns2Write
    mcols(fds, type="psi5") <- cbind(
            mcols(fds, type="psi5")[,col2Take], newAnno[,..columns2Write])
    
    # return annotated object
    fds
}


