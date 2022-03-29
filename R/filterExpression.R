#' @title Filtering FraserDataSets
#' 
#' @description This method can be used to filter out introns that are not 
#' reliably detected and to remove introns with no variablity between samples.
#' 
#' @inheritParams countRNA
#' @param object A \code{\link{FraserDataSet}} object
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
#' @return A FraserDataSet with information about which junctions passed the
#' filters. If \code{filter=TRUE}, the filtered FraserDataSet is returned.
#' 
#' @examples
#' fds <- createTestFraserDataSet()
#' fds <- filterExpressionAndVariability(fds, minDeltaPsi=0.1, filter=FALSE)
#' mcols(fds, type="psi5")[, c(
#'         "maxCount", "passedExpression", "maxDPsi3", "passedVariability")]
#' 
#' plotFilterExpression(fds)
#' plotFilterVariability(fds)     
#' 
#' @name filtering
#' @rdname filtering
#' 
NULL

#' @describeIn filtering This functions filters out both introns with low 
#' read support and introns that are not variable across samples. 
#' @export
filterExpressionAndVariability <- function(object, minExpressionInOneSample=20, 
                    quantile=0.05, quantileMinExpression=1, minDeltaPsi=0,
                    filter=TRUE, 
                    delayed=ifelse(ncol(object) <= 300, FALSE, TRUE),
                    BPPARAM=bpparam()){
    # filter introns with low read support and corresponding splice sites
    object <- filterExpression(object, 
                    minExpressionInOneSample=minExpressionInOneSample, 
                    quantile=quantile, 
                    quantileMinExpression=quantileMinExpression, 
                    filter=filter, delayed=delayed,
                    BPPARAM=BPPARAM)
    
    # filter introns that are not variable across samples
    object <- filterVariability(object, minDeltaPsi=minDeltaPsi, filter=filter, 
                    delayed=delayed, BPPARAM=BPPARAM)
    
    # return fds
    message(date(), ": Filtering done!")
    return(object)
    
}

filterExpression.FRASER <- function(object, minExpressionInOneSample=20,
                    quantile=0.05, quantileMinExpression=1, filter=TRUE, 
                    delayed=ifelse(ncol(object) <= 300, FALSE, TRUE),
                    BPPARAM=bpparam()){

    stopifnot(is(object, "FraserDataSet"))
    
    message(date(), ": Filtering out introns with low read support ...")
    
    # extract counts
    cts  <- K(object, type="j")
    ctsN5 <- N(object, type="psi5")
    ctsN3 <- N(object, type="psi3")
    
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
        mcols(object, type="j")[n] <- cutoffs[[n]]
    }
    
    mcols(object, type="j")[['passedExpression']] <-
            cutoffs$maxCount >= minExpressionInOneSample &
            (cutoffs$quantileValue5 >= quantileMinExpression |
                cutoffs$quantileValue3 >= quantileMinExpression) 
    if("passedVariability" %in% colnames(mcols(object, type="j"))){
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedExpression']] & 
            mcols(object, type="j")[['passedVariability']]
    } else{
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedExpression']]
    }
    
    # filter if requested
    if(isTRUE(filter)){
        object <- applyExpressionFilters(object, minExpressionInOneSample, 
                quantileMinExpression)
    }

    validObject(object)
    return(object)
}

#' @describeIn filtering This function filters out introns and corresponding 
#' splice sites that have low read support in all samples.
#' @export
setMethod("filterExpression", signature="FraserDataSet",
        filterExpression.FRASER)

#' @describeIn filtering This function filters out introns and corresponding 
#' splice sites which do not show variablity across samples.
#' @export
filterVariability <- function(object, minDeltaPsi=0, filter=TRUE, 
                    delayed=ifelse(ncol(object) <= 300, FALSE, TRUE),
                    BPPARAM=bpparam()){
    
    message(date(), ": Filtering out non-variable introns ...")
    
    # extract counts
    cts    <- K(object, type="j")
    ctsSE  <- K(object, type="ss")
    ctsN5  <- N(object, type="psi5")
    ctsN3  <- N(object, type="psi3")
    ctsNSE <- N(object, type="theta")
    
    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        ctsN5 <- as.matrix(ctsN5)
        ctsN3 <- as.matrix(ctsN3)
        ctsSE <- as.matrix(ctsSE)
        ctsNSE <- as.matrix(ctsNSE)
    }
    
    # cutoff functions
    f1 <- function(cts, ctsN3, ...) {
        psi <- cts/ctsN3
        rowMaxs(abs(psi - rowMeans2(psi, na.rm=TRUE)), na.rm=TRUE) }
    f2 <- function(cts, ctsN5, ...) {
        psi <- cts/ctsN5
        rowMaxs(abs(psi - rowMeans2(psi, na.rm=TRUE)), na.rm=TRUE) }
    f3 <- function(ctsSE, ctsNSE, ...) {
        theta <- ctsSE/ctsNSE
        dTheta <- rowMaxs(abs(theta - rowMeans2(theta, na.rm=TRUE)), 
                            na.rm=TRUE) }
        
    
    funs <- c(maxDPsi3=f1, maxDPsi5=f2, maxDTheta=f3)
    
    # run it in parallel
    cutoffs <- bplapply(funs, function(f, ...) f(...), BPPARAM=BPPARAM,
                        cts=cts, ctsN3=ctsN3, ctsN5=ctsN5, 
                        ctsSE=ctsSE, ctsNSE=ctsNSE)
    
    # add annotation to object
    for(n in names(cutoffs)){
        if(n == "maxDTheta"){
            mcols(object, type="ss")[n] <- cutoffs[[n]]
        } else{ 
            mcols(object, type="j")[n] <- cutoffs[[n]]
        }
    }
    
    # add annotation of theta on splice sites of introns to mcols
    intron_dt <- as.data.table(rowRanges(object, type="j"))
    ss_dt <- as.data.table(rowRanges(object, type="ss"))
    mcols(object, type="j")["maxDThetaDonor"] <- 
        merge(intron_dt, ss_dt, by.x="startID", by.y="spliceSiteID", 
                all.x=TRUE, sort=FALSE)[,maxDTheta]
    mcols(object, type="j")["maxDThetaAcceptor"] <- 
        merge(intron_dt, ss_dt, by.x="endID", by.y="spliceSiteID", 
                all.x=TRUE, sort=FALSE)[,maxDTheta]

    # check which introns pass the filter
    mcols(object, type="j")[['passedVariability']] <- pmax(na.rm=TRUE,
            cutoffs$maxDPsi3, 
            cutoffs$maxDPsi5, 
            mcols(object, type="j")$maxDThetaDonor, 
            mcols(object, type="j")$maxDThetaAcceptor,
            0) >= minDeltaPsi
    if("passedExpression" %in% colnames(mcols(object, type="j"))){
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedExpression']] & 
            mcols(object, type="j")[['passedVariability']]
    } else{
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedVariability']]
    }
    
    # filter if requested
    if(isTRUE(filter)){
        object <- applyVariabilityFilters(object, minDeltaPsi)
    }
    
    validObject(object)
    return(object)
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
    numFilt <- sum(mcols(fds, type="j")[['passedExpression']])
    message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                    ". This is ", signif(numFilt/length(fds)*100, 3),
                    "% of the junctions"))
    fds <- fds[mcols(fds, type="j")[['passedExpression']], by="psi5"]
    
    return(fds)
    
}

#' Applies previously calculated variablilty filters
#' @noRd
applyVariabilityFilters <- function(fds, minDeltaPsi){
    
    #
    passedVariability <-  mcols(fds, type="j")[['passedVariability']]
    # maxDPsi3 <- mcols(fds, type="j")[['maxDPsi3']]
    # maxDPsi5 <- mcols(fds, type="j")[['maxDPsi5']]
    # maxDThetaDonor    <- mcols(fds, type="j")[['maxDThetaDonor']]
    # maxDThetaAcceptor <- mcols(fds, type="j")[['maxDThetaAcceptor']]
    
    # store information of non-variable junctions
    filtered <- !passedVariability
    # filtered <- (pmax(maxDPsi3, maxDPsi5, maxDThetaDonor, maxDThetaAcceptor) 
    #                 < minDeltaPsi)
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
                        mcols(nV_stored)$maxDPsi3 < minDeltaPsi &
                        mcols(nV_stored)$maxDThetaDonor < minDeltaPsi &
                        mcols(nV_stored)$maxDThetaAcceptor < minDeltaPsi
            
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
    numFilt <- sum(passedVariability)
    message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                    ". This is ", signif(numFilt/length(fds)*100, 3),
                    "% of the junctions"))
    fds <- fds[mcols(fds, type="j")[['passedVariability']], by="psi5"]
    return(fds)
        
}

#' @describeIn filtering This functions filters out both introns with low 
#' read support and introns that are not variable across samples. 
#' @export
filterExpressionAndVariability_jaccard <- function(object, 
                        minExpressionInOneSample=20, 
                        quantile=0.95, quantileMinExpression=1, minDelta=0.05,
                        filter=TRUE, 
                        delayed=ifelse(ncol(object) <= 300, FALSE, TRUE),
                        BPPARAM=bpparam()){
    # filter introns with low read support and corresponding splice sites
    object <- filterExpression_jaccard(object, 
                        minExpressionInOneSample=minExpressionInOneSample, 
                        quantile=quantile, 
                        quantileMinExpression=quantileMinExpression, 
                        filter=filter, delayed=delayed,
                        BPPARAM=BPPARAM)
    
    # filter introns that are not variable across samples
    object <- filterVariability_jaccard(object, minDelta=minDelta, 
                        filter=filter, 
                        delayed=delayed, BPPARAM=BPPARAM)
    
    # return fds
    message(date(), ": Filtering done!")
    return(object)
    
}

#' @describeIn filtering This function filters out introns and corresponding 
#' splice sites which are expressed at very low levels across samples.
#' @export
filterExpression_jaccard <- function(object, minExpressionInOneSample=20,
                            quantile=0.95, quantileMinExpression=1, filter=TRUE, 
                            delayed=ifelse(ncol(object) <= 300, FALSE, TRUE),
                            BPPARAM=bpparam()){
    
    stopifnot(is(object, "FraserDataSet"))
    
    message(date(), ": Filtering out introns with low read support ...")
    
    # extract counts
    cts  <- K(object, type="j")
    ctsN <- N(object, type="jaccard")
    
    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        ctsN <- as.matrix(ctsN)
    }
    
    # cutoff functions
    f1 <- function(cts, ...){
        rowMaxs(cts) }
    f2 <- function(cts, ctsN, quantile, ...){
        rowQuantiles(ctsN, probs=quantile, drop=FALSE)[,1] }
    
    funs <- c(maxCount=f1, quantileValueN=f2)
    
    # run it in parallel
    cutoffs <- bplapply(funs, function(f, ...) f(...), BPPARAM=BPPARAM,
                        cts=cts, ctsN=ctsN, quantile=quantile)
    
    # add annotation to object
    for(n in names(cutoffs)){
        mcols(object, type="j")[n] <- cutoffs[[n]]
    }
    
    mcols(object, type="j")[['passedExpression']] <-
        cutoffs$maxCount >= minExpressionInOneSample &
        cutoffs$quantileValueN >= quantileMinExpression  
    if("passedVariability" %in% colnames(mcols(object, type="j"))){
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedExpression']] & 
            mcols(object, type="j")[['passedVariability']]
    } else{
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedExpression']]
    }
    
    # filter if requested
    if(isTRUE(filter)){
        object <- applyExpressionFilters_jaccard(object, 
                                                minExpressionInOneSample, 
                                                quantileMinExpression)
    }
    
    validObject(object)
    return(object)
}

#' @describeIn filtering This function filters out introns and corresponding 
#' splice sites which do not show variablity across samples.
#' @export
filterVariability_jaccard <- function(object, minDelta=0, filter=TRUE, 
                              delayed=ifelse(ncol(object) <= 300, FALSE, TRUE),
                              BPPARAM=bpparam()){
    
    message(date(), ": Filtering out non-variable introns ...")
    
    # extract counts
    cts    <- K(object, type="j")
    ctsN  <- N(object, type="jaccard")
    
    if(isFALSE(delayed)){
        cts <- as.matrix(cts)
        ctsN <- as.matrix(ctsN)
    }
    
    # cutoff functions
    f1 <- function(cts, ctsN, ...) {
        jaccard <- cts/ctsN
        rowMaxs(abs(jaccard - rowMeans2(jaccard, na.rm=TRUE)), 
                na.rm=TRUE) }
    
    funs <- c(maxDJaccard=f1)
    
    # run it in parallel
    cutoffs <- bplapply(funs, function(f, ...) f(...), BPPARAM=BPPARAM,
                        cts=cts, ctsN=ctsN)
    
    # add annotation to object
    for(n in names(cutoffs)){
            mcols(object, type="j")[n] <- cutoffs[[n]]
    }
    
    # add annotation of theta on splice sites of introns to mcols
    intron_dt <- as.data.table(rowRanges(object, type="j"))
    
    # check which introns pass the filter
    mcols(object, type="j")[['passedVariability']] <- pmax(na.rm=TRUE,
                                                           cutoffs$maxDJaccard, 
                                                           0) >= minDelta
    if("passedExpression" %in% colnames(mcols(object, type="j"))){
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedExpression']] & 
            mcols(object, type="j")[['passedVariability']]
    } else{
        mcols(object, type="j")[['passed']] <-  
            mcols(object, type="j")[['passedVariability']]
    }
    
    # filter if requested
    if(isTRUE(filter)){
        object <- applyVariabilityFilters(object, minDelta)
    }
    
    validObject(object)
    return(object)
}

#' Applies previously calculated filters for expression filters
#' @noRd
applyExpressionFilters_jaccard <- function(fds, minExpressionInOneSample, 
                                   quantileMinExpression){
    
    maxCount       <- mcols(fds, type="j")[['maxCount']]
    quantileValueN <- mcols(fds, type="j")[['quantileValueN']]
    
    # report rare junctions that passed minExpression filter but not 
    # quantileFilter as SE obj
    junctionsToReport <- maxCount >= minExpressionInOneSample & 
        !(quantileValueN >= quantileMinExpression) 
    outputDir <- file.path(workingDir(fds), "savedObjects", nameNoSpace(fds))
    
    if(any(junctionsToReport)){
        # get SE object of junctions to report
        rareJunctions <- asSE(fds[junctionsToReport, by="j"])
        for(aname in assayNames(rareJunctions)){
            if(!(aname %in% c("rawCountsJ", "rawOtherCounts_psi5", 
                              "rawOtherCounts_psi3", "psi5", "psi3", 
                              "delta_psi5", "delta_psi3", "jaccard",
                              "rawOtherCounts_intron_jaccard"))){
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
                !(mcols(rJ_stored)$quantileValueN >= quantileMinExpression) 
            
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
    numFilt <- sum(mcols(fds, type="j")[['passedExpression']])
    message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                   ". This is ", signif(numFilt/length(fds)*100, 3),
                   "% of the junctions"))
    fds <- fds[mcols(fds, type="j")[['passedExpression']], by="psi5"]
    
    return(fds)
    
}


#' Applies previously calculated variablilty filters
#' @noRd
applyVariabilityFilters <- function(fds, minDelta){
    
    #
    passedVariability <-  mcols(fds, type="j")[['passedVariability']]
    
    # store information of non-variable junctions
    filtered <- !passedVariability

    outputDir <- file.path(workingDir(fds), "savedObjects", nameNoSpace(fds))
    if(any(filtered)){
        # get SE object of junctions to report
        nonVariableJunctions <- asSE(fds[filtered, by="j"])
        for(aname in assayNames(nonVariableJunctions)){
            if(!(aname %in% c("rawCountsJ", "rawOtherCounts_psi5", 
                              "rawOtherCounts_psi3", "psi5", "psi3", 
                              "delta_psi5", "delta_psi3", "jaccard",
                              "rawOtherCounts_intron_jaccard"))){
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
            toReport <- mcols(nV_stored)$maxDJaccard < minDelta
            
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
    numFilt <- sum(passedVariability)
    message(paste0("Keeping ", numFilt, " junctions out of ", length(fds),
                   ". This is ", signif(numFilt/length(fds)*100, 3),
                   "% of the junctions"))
    fds <- fds[mcols(fds, type="j")[['passedVariability']], by="psi5"]
    return(fds)
    
}
