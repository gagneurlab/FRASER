#'
#' Annotates the given FraserDataSet with the HGNC symbol with biomaRt
#' 
#' @param fds FraserDataSet
#' @param feature Defines which feature (default is HGNC symbol) should be 
#' annotated.
#' @param featureName Name of the feature in the FraserDataSet mcols.
#' @param biotype The biotype.
#' @param ensembl The ensembl that should be used. If NULL, the default one is 
#' used (hsapiens_gene_ensembl, GRCh37).
#' @param GRCh GRCh version to connect to. If this is NULL, then the current 
#' GRCh38 is used. Otherwise, this can only be 37 (default) at the moment 
#' (see \code{\link[biomaRt]{useEnsembl}}).
#' @param txdb A \code{TxDb} object. If this is NULL, then the default 
#' one is used, currently this is \code{TxDb.Hsapiens.UCSC.hg19.knownGene}.
#' @param orgDb An \code{orgDb} object If this is NULL, then the 
#' default one is used, currently this is \code{org.Hs.eg.db}.
#' @param keytype The type of gene IDs in the \code{TxDb} object (see 
#' AnnotationDbi::keytypes(orgDb) for a list of available ID types).
#' 
#' @return FraserDataSet
#' 
#' @examples
#'
#' fds <- createTestFraserDataSet()
#' 
#' ### Two ways to annotage ranges with gene names: 
#' # either using biomart:
#' fds <- annotateRanges(fds, GRCh=NULL)
#' rowRanges(fds, type="psi5")[,"hgnc_symbol"]
#'  
#' # or with a TxDb object
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
    stopifnot(is(fds, "FraserDataSet"))
    if(length(fds) == 0) return(fds)
    
    # useEnsembl only understands GRCh=37 or GRCh=NULL (uses 38 then)
    if(GRCh == 38){
        GRCh <- NULL
    }

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
        mcols(fds, type=i)[[featureName]] <- annos[["feature"]]
        mcols(fds, type=i)[[paste0("other_", featureName)]] <- 
            annos[["other_features"]]
    }

    return(fds)
}

#' @rdname annotateRanges
#' @export
annotateRangesWithTxDb <- function(fds, feature="SYMBOL", 
                                    featureName="hgnc_symbol",
                                    keytype="ENTREZID",
                                    txdb=NULL, orgDb=NULL){
    
    # check input
    stopifnot(is(fds, "FraserDataSet"))
    if(length(fds) == 0) return(fds)
    
    if(is.null(txdb)){
        if(requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")){
            txdb <- 
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene 
        } else{
            stop("Please provide a TxDb object as input.")
        }
        
    }
    if(is.null(orgDb)){
        if(requireNamespace("org.Hs.eg.db")){
            orgDb <- org.Hs.eg.db::org.Hs.eg.db
        } else{
            stop("Please provide an OrgDb object to extract gene symbols")
        }
    }
    
    for(i in c("psi3", "psiSite")){
        # get GRanges object with the split reads which should be annotated
        gr <- rowRanges(fds, type=i)
        
        # get the annotation to compare to
        anno <- genes(txdb)
        tmp  <- select(orgDb, keys=mcols(anno)[,"gene_id"], columns=feature, 
                    keytype=keytype)
        tmp  <- as.data.table(tmp)
        tmp[, uniqueID := .GRP, by=keytype]
        anno <- anno[tmp[,uniqueID]]
        mcols(anno)[[featureName]] <- tmp[,get(feature)]
        anno <- anno[!is.na(mcols(anno)[,featureName]),]
        anno <- anno[mcols(anno)[,featureName] != "",]
        if(any(strand(gr) == "*")){
            strand(anno) <- "*"
        }
        
        # retrieve the feature of interest for the split reads
        annos <- getAnnotationFeature(data=gr, featureName, anno)
        mcols(fds, type=i)[[featureName]] <- annos[["feature"]]
        mcols(fds, type=i)[[paste0("other_", featureName)]] <- 
            annos[["other_features"]]
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
    # featureDT <- unique(featureDT, by="from")
    featureDT <- featureDT[,
            list(first_feature=unique(feature)[1], 
                other_features=paste(unique(feature)[-1], collapse = ";")),
            by="from"
    ]

    return(list(feature=featureDT[order(from),first_feature],
                other_features=featureDT[order(from),other_features]))
}


#'
#' Find annotated junctions
#' 
#' Annotate the object with a given annotation. 
#' 
#' @return FraserDataSet
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
    
    if(isFALSE(as.logical(stranded))){
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


