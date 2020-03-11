#' @title Filtering FraseRDataSets
#' 
#' @description This method can be used to filter out introns that are not 
#' reliably detected and to remove introns with no variablity between samples.
#' 
#' @inheritParams countRNA
#' @param minExpressionInOneSample The minimal read count in at least one 
#' sample that is required for an intron to pass the filter.
#' @param quantile Defines which quantile should be considered for the filter.
#' @param quantileMinExpression The minimum read count an intron needs to have 
#' at the specified quantile to pass the filter.
#' @param minDeltaPsi Only introns for which the maximal difference in the psi 
#' value of a sample to the mean psi of the intron is larger than this value 
#' pass the filter. 
#' @param filter If TRUE, a subsetted fds containing only the introns that 
#' passed all filters is returned. If FALSE, no subsetting is done and the 
#' information of whether an intron passed the filters is only stored in the 
#' mcols.
#' @param delayed If FALSE, count matrices will be loaded into memory, 
#' otherwise the function works on the delayedMatrix representations. The 
#' default value depends on the number of samples in the fds-object. 
#'
#' @return A FraseRDataSet with information about which junctions passed the
#' filters. If \code{filter=TRUE}, the filtered FraseRDataSet is returned.
#' 
#' @examples
#' fds <- makeExampleFraseRDataSet()
#' fds <- calculatePSIValues(fds)
#' fds <- filterExpressionAndVariability(fds)
#'
#' @name filtering
#' @rdname filtering
#' 
NULL

#' @describeIn filtering This functions filters out both introns with low 
#' read support and introns that are not variable across samples. 
#' @export
filterExpressionAndVariability <- function(fds, minExpressionInOneSample=20, 
                                quantile=0.05, quantileMinExpression=1, 
                                minDeltaPsi=0, filter=TRUE, 
                                delayed=ifelse(ncol(fds) <= 300, FALSE, TRUE),
                                BPPARAM=bpparam()){
    # filter introns with low read support and corresponding splice sites
    fds <- filterExpression(fds, 
                    minExpressionInOneSample=minExpressionInOneSample, 
                    quantile=quantile, 
                    quantileMinExpression=quantileMinExpression, 
                    filter=filter, delayed=delayed,
                    BPPARAM=BPPARAM)
    
    # filter introns that are not variable across samples
    fds <- filterVariability(fds, minDeltaPsi=minDeltaPsi, filter=filter, 
                    delayed=delayed, BPPARAM=BPPARAM)
    
    # return fds
    message(date(), ": Filtering done!")
    return(fds)
    
}

#' @describeIn filtering This function filters out introns and corresponding 
#' splice sites that have low read support in all samples.
#' @export
filterExpression <- function(fds, minExpressionInOneSample=20, quantile=0.05,
                    quantileMinExpression=1, filter=TRUE, 
                    delayed=ifelse(ncol(fds) <= 300, FALSE, TRUE),
                    BPPARAM=bpparam()){

    stopifnot(is(fds, "FraseRDataSet"))
    
    message(date(), ": Filtering out introns with low read support ...")
    
    # extract counts
    cts  <- K(fds, type="j")
    ctsN5 <- N(fds, type="psi5")
    ctsN3 <- N(fds, type="psi3")
    
    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        ctsN5 <- as.matrix(ctsN5)
        ctsN3 <- as.matrix(ctsN3)
    }

    # cutoff functions
    f1 <- function(cts, ...){
            rowMaxs(cts) }
    f2 <- function(cts, ctsN5, quantile, ...){
            rowQuantiles(ctsN5, probs=quantile, drop=FALSE)[,1] }
    f3 <- function(cts, ctsN3, quantile, ...) {
            rowQuantiles(ctsN3, probs=quantile, drop=FALSE)[,1] }

    funs <- c(maxCount=f1, quantileValue5=f2, quantileValue3=f3)

    # run it in parallel
    cutoffs <- bplapply(funs, function(f, ...) f(...), BPPARAM=BPPARAM,
            cts=cts, ctsN3=ctsN3, ctsN5=ctsN5, quantile=quantile)

    # add annotation to object
    for(n in names(cutoffs)){
        mcols(fds, type="j")[n] <- cutoffs[[n]]
    }
    
    mcols(fds, type="j")[['passed']] <-
            cutoffs$maxCount >= minExpressionInOneSample &
            (cutoffs$quantileValue5 >= quantileMinExpression |
                cutoffs$quantileValue3 >= quantileMinExpression) 
    
    # filter if requested
    if(isTRUE(filter)){
        fds <- applyExpressionFilters(fds, minExpressionInOneSample, 
                                quantileMinExpression)
    }

    validObject(fds)
    return(fds)
}

#' @describeIn filtering This function filters out introns and corresponding 
#' splice sites which do not show variablity across samples.
#' @export
filterVariability <- function(fds, minDeltaPsi=0, filter=TRUE, 
                                delayed=ifelse(ncol(fds) <= 300, FALSE, TRUE),
                                BPPARAM=bpparam()){
    
    message(date(), ": Filtering out non-variable introns ...")
    
    # extract counts
    cts  <- K(fds, type="j")
    ctsN5 <- N(fds, type="psi5")
    ctsN3 <- N(fds, type="psi3")
    
    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        ctsN5 <- as.matrix(ctsN5)
        ctsN3 <- as.matrix(ctsN3)
    }
    
    # cutoff functions
    f1 <- function(cts, ctsN3, ...) {
        psi <- cts/ctsN3
        rowMaxs(abs(psi - rowMeans2(psi, na.rm=TRUE)), na.rm=TRUE) }
    f2 <- function(cts, ctsN5, ...) {
        psi <- cts/ctsN5
        rowMaxs(abs(psi - rowMeans2(psi, na.rm=TRUE)), na.rm=TRUE) }
    
    funs <- c(maxDPsi3=f1, maxDPsi5=f2)
    
    # run it in parallel
    cutoffs <- bplapply(funs, function(f, ...) f(...), BPPARAM=BPPARAM,
                        cts=cts, ctsN3=ctsN3, ctsN5=ctsN5)
    
    # add annotation to object
    for(n in names(cutoffs)){
        mcols(fds, type="j")[n] <- cutoffs[[n]]
    }
    mcols(fds, type="j")[['passed']] <-
        (cutoffs$maxDPsi3     >= minDeltaPsi |
            cutoffs$maxDPsi5 >= minDeltaPsi)
    
    # filter if requested
    if(isTRUE(filter)){
        fds <- applyVariabilityFilters(fds, minDeltaPsi)
    }
    
    validObject(fds)
    return(fds)
}

#' Applies previously calculated filters for expression filters
#' @noRd
applyExpressionFilters <- function(fds, minExpressionInOneSample, 
                                    quantileMinExpression){
    
    maxCount       <- mcols(fds, type="j")[['maxCount']]
    quantileValue5 <- mcols(fds, type="j")[['quantileValue5']]
    quantileValue3 <- mcols(fds, type="j")[['quantileValue3']]
    
    # report rare junctions that passed minExpression filter but not 
    # quantileFilter as SE obj
    junctionsToReport <- maxCount >= minExpressionInOneSample & 
                        !(quantileValue5 >= quantileMinExpression |
                                quantileValue3 >= quantileMinExpression) 
    outputDir <- file.path(workingDir(fds), "savedObjects", nameNoSpace(fds))
    
    if(any(junctionsToReport)){
        # get SE object of junctions to report
        rareJunctions <- asSE(fds[junctionsToReport, by="j"])
        for(aname in assayNames(rareJunctions)){
            if(!(aname %in% c("rawCountsJ", "rawOtherCounts_psi5", 
                                "rawOtherCounts_psi3", "psi5", "psi3", 
                                "delta_psi5", "delta_psi3"))){
                assay(rareJunctions, aname) <- NULL
            }
        }
        rareJunctions <- saveHDF5SummarizedExperiment(rareJunctions, 
                                            dir=file.path(tempdir(), "tmp_rJ"), 
                                            replace=TRUE)
        
        # check if folder already exists from previous filtering
        rareJctsDir <- file.path(outputDir, "rareJunctions")
        if(dir.exists(rareJctsDir)){
            warning("Filtering has already been applied previously. Introns ", 
                    "that were already filtered out but should be kept now ",
                    "cannot be restored.")
            rJ_stored <- loadHDF5SummarizedExperiment(dir=rareJctsDir)
            toReport <- mcols(rJ_stored)$maxCount >= minExpressionInOneSample & 
                !(mcols(rJ_stored)$quantileValue5 >= quantileMinExpression |
                    mcols(rJ_stored)$quantileValue3 >= quantileMinExpression) 
            
            rJ_tmp <- rbind(rJ_stored[toReport,], rareJunctions)
            
            for(aname in assayNames(rJ_tmp)){
                assay(rJ_tmp, aname) <- 
                    rbind(as.matrix(assay(rareJunctions, aname)), 
                            as.matrix(assay(rJ_stored[toReport,], aname)) )
            }
            rareJunctions <- rJ_tmp
            rm(rJ_tmp)
        } 
        
        rareJunctions <- saveHDF5SummarizedExperiment(rareJunctions, 
                                                dir=rareJctsDir, replace=TRUE)
    }
    
    # apply filter
    numFilt <- sum(mcols(fds, type="j")[['passed']])
    message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                    ". This is ", signif(numFilt/length(fds)*100, 3),
                    "% of the junctions"))
    fds <- fds[mcols(fds, type="j")[['passed']], by="psi5"]
    
    return(fds)
    
}

#' Applies previously calculated variablilty filters
#' @noRd
applyVariabilityFilters <- function(fds, minDeltaPsi){
    
    #
    maxDPsi3 <- mcols(fds, type="j")[['maxDPsi3']]
    maxDPsi5 <- mcols(fds, type="j")[['maxDPsi5']]
    
    # store information of non-variable junctions
    filtered <- !(maxDPsi3 >= minDeltaPsi | maxDPsi5 >= minDeltaPsi)
    outputDir <- file.path(workingDir(fds), "savedObjects", nameNoSpace(fds))
    if(any(filtered)){
        # get SE object of junctions to report
        nonVariableJunctions <- asSE(fds[filtered, by="j"])
        for(aname in assayNames(nonVariableJunctions)){
            if(!(aname %in% c("rawCountsJ", "rawOtherCounts_psi5", 
                                "rawOtherCounts_psi3", "psi5", "psi3", 
                                "delta_psi5", "delta_psi3"))){
                assay(nonVariableJunctions, aname) <- NULL
            }
        }
        nonVariableJunctions <- saveHDF5SummarizedExperiment(replace=TRUE,
                                        nonVariableJunctions,
                                        dir=file.path(tempdir(), "tmp_nvJ"))
        
        # check if folder already exists from previous filtering
        nonVarJctsDir <- file.path(outputDir, "nonVariableJunctions")
        if(dir.exists(nonVarJctsDir)){
            warning("Filtering has already been applied previously. Introns ", 
                    "that were already filtered out but should be kept now ",
                    "cannot be restored.")
            nV_stored <- loadHDF5SummarizedExperiment(dir=nonVarJctsDir) 
            toReport <- mcols(nV_stored)$maxDPsi5 < minDeltaPsi &
                        mcols(nV_stored)$maxDPsi3 < minDeltaPsi
            
            nVJunctions <- rbind(nonVariableJunctions, nV_stored[toReport,])
            for(aname in assayNames(nVJunctions)){
                assay(nVJunctions, aname) <- 
                rbind(as.matrix(assay(nonVariableJunctions, aname)), 
                        as.matrix(assay(nV_stored[toReport,], aname)) )
            }
            nonVariableJunctions <- nVJunctions
            rm(nVJunctions)
        } 
        
        nonVariableJunctions <- saveHDF5SummarizedExperiment(dir=nonVarJctsDir,
                                        x=nonVariableJunctions, replace=TRUE)
        
    }
    
    # apply filtering
    numFilt <- sum(mcols(fds, type="j")[['passed']])
    message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                    ". This is ", signif(numFilt/length(fds)*100, 3),
                    "% of the junctions"))
    fds <- fds[mcols(fds, type="j")[['passed']], by="psi5"]
    return(fds)
        
}
