########
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for annotating the ranges with
## biomaRt from ENSEMBL
##

#'
#' Annotates the given FraseRDataSet with the HGNC symbol with biomaRt
#' 
#' @param fds FraseRDataSet
#' @param feature Defines which feature (default is HGNC symbol) should be 
#' annotated.
#' @param featureName Name of the feature in the FraseRDataSet mcols.
#' @param biotype The biotype.
#' @param ensembl The ensembl that should be used. If NULL, the default one is 
#' used (hsapiens_gene_ensembl, GRCh37).
#' @param GRCh GRCh version to connect to. If this is NULL, then the current 
#' GRCh38 is used. Otherwise, this can only be 37 (default) at the moment 
#' (see \code{\link{useEnsebml}}).
#' @param txdb A TxDb object (default: TxDb.Hsapiens.UCSC.hg19.knownGene).
#' @param orgDb An orgDb object (default: org.Hs.eg.db).
#' 
#' @return FraseRDataSet
#' 
#' @examples
#'
#' fds <- makeExampleFraseRDataSet()
#' fds <- annotateRanges(fds, GRCh=NULL)
#' 
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' require(org.Hs.eg.db)
#' orgDb <- org.Hs.eg.db
#' fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
#'
#' rowRanges(fds, type="psi5")[,"hgnc_symbol"]
#'
#' @rdname annotateRanges
#' @export
annotateRanges <- function(fds, feature="hgnc_symbol", featureName=feature,
            biotype=list("protein_coding"), ensembl=NULL, GRCh=37){

    # check input
    stopifnot(is(fds, "FraseRDataSet"))
    if(length(fds) == 0) return(fds)

    if(is.null(ensembl)){
        tryCatch({
            ensemblOutput <- capture.output(ensembl <- useEnsembl(
                    biomart="ensembl", dataset="hsapiens_gene_ensembl", 
                    GRCh=GRCh
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

#' @rdname annotateRanges
#' @export
annotateRangesWithTxDb <- function(fds, feature="SYMBOL", 
                                    featureName="hgnc_symbol",
                                    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                    orgDb=org.Hs.eg.db){
 
    # check input
    stopifnot(is(fds, "FraseRDataSet"))
    if(length(fds) == 0) return(fds)
    
    for(i in c("psi3", "psiSite")){
        # get GRanges object with the split reads which should be annotated
        gr <- rowRanges(fds, type=i)
        
        # get the annotation to compare to
        anno <- genes(txdb)
        mcols(anno)[[featureName]] <- 
            select(orgDb, keys=mcols(anno)[,"gene_id"], columns=feature, 
                   keytype="ENTREZID")[,feature]
        anno <- anno[!is.na(mcols(anno)[,featureName]),]
        anno <- anno[mcols(anno)[,featureName] != "",]
        if(any(strand(gr) == "*")){
            strand(anno) <- "*"
        }
        
        # retrieve the feature of interest for the split reads
        mcols(fds, type=i)[[featureName]] <- 
            getAnnotationFeature(gr, featureName, anno)
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
#' @examples 
#'     TODO <- 1
#' @noRd
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


