##
## @author Christian Mertes \email{mertes@@in.tum.de}
## 
## This file contains all functions for calculating the statistics
## First it starts with calculating the Z-score for each
## site and then the p-values are calculated dependend on the 
## given method in the setting file 
##

#'
#' calculate the zscore for each psi value
#' 
#' @export
calculateZScores <- function(dataset){
    
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")
    
    message(date(), ": Calculate the PSI3 values ...")
    dataset <- .calculateZScorePerDataSet(dataset, "splitReads", "psi3")
    message(date(), ": Calculate the PSI5 values ...")
    dataset <- .calculateZScorePerDataSet(dataset, "splitReads", "psi5")
    message(date(), ": Calculate the sitePSI values ...")
    dataset <- .calculateZScorePerDataSet(dataset, "nonSplicedReads", "sitePSI")
    
    return(dataset)
}

#'
#' calculates the zscore for a given data type and a given psi type
#' and adds it directly to the dataset itself
#'
#' @noRd
.calculateZScorePerDataSet <- function(dataset, readType, psiType){
    
    # data to work with
    seCounts <- slot(dataset, readType)
  
    # get raw data and replace NA's with zeros
    psiVal <- .getAssayAsDataTable(seCounts, psiType)
    
    # z = ( x - mean ) / sd
    rowmean <- rowMeans(psiVal, na.rm = TRUE)
    rowsd   <- apply(psiVal, 1, sd, na.rm = TRUE)
    zscores <- (psiVal - rowmean) / rowsd
    
    # add it to the FraseR object
    assayName <- paste0("zscore_", psiType)
    zscores   <- .asDataFrame(zscores, dataset@settings@sampleData[,sampleID])
    assays(slot(dataset, readType))[[assayName]] <- zscores
          
    return(dataset)  
}

#'
#' calculates the P-Value for the given FraseR dataset object
#' The P-Value calculation is based on the given method in the
#' FraseRSettings object
#'
#' @export
calculatePValues <- function(dataset, internBPPARAM=SerialParam()){
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")
    
    # check which method we should use
    method <- dataset@settings@method
    if(method == "Fisher"){
        return(.testPsiWithFisher(dataset, internBPPARAM))
    }
    
    if(method == "DESeq2"){
        stop("This method is not yet implemented.")
    }
    
    if(method == "Martin"){
        stop("This method is not yet implemented.")
    }
    
    stop("The provided method is not present for this package.",
            "Please set the method to one of the following: Fisher, DESeq2, Martin"
    )
}

#'
#'
#' @noRd
.testPsiWithFisher <- function(dataset, internBPPARAM){
    
    # test all 3 different types
    assays(dataset@splitReads)$pvalue_psi3 <- 
            .testPsiWithFisherPerType(dataset, "splitReads", "psi3", internBPPARAM)
    assays(dataset@splitReads)$pvalue_psi5 <- 
            .testPsiWithFisherPerType(dataset, "splitReads", "psi5", internBPPARAM)
    assays(dataset@nonSplicedReads)$pvalue_sitePSI <- 
            .testPsiWithFisherPerType(dataset, "nonSplicedReads", "sitePSI", internBPPARAM)
    
    # return the new datasets
    return(dataset)
}

#'
#'
#' @noRd
.testPsiWithFisherPerType <- function(dataset, readType, psiType, internBPPARAM){
    # go over each group but no NA's
    group   <- sampleGroup(dataset@settings)
    
    # reads to test for abberent splicing (eg: nonSplicedReads)
    rawCounts <- .getAssayAsDataTable(slot(dataset, readType), "rawCounts")
    
    # other reads (eg: splitReads)
    rawOtherCounts <- .getAssayAsDataTable(slot(dataset, readType), paste0("rawOtherCounts_", psiType))
    
    pvalues <- bplapply(unique(na.omit(group)), dataset=dataset, 
                        rawCounts=rawCounts, rawOtherCounts=rawOtherCounts,
                        BPPARAM=dataset@settings@parallel,
                        internBPPARAM=internBPPARAM,
                        FUN=.testPsiWithFisherPerGroup
    )
    names(pvalues) <- as.character(unique(na.omit(group)))
    pvalues_full <- pvalues[as.character(group)]
    
    # add NA's to the non tested ones
    pvalues_full[is.na(group)] <- list(rep(as.numeric(NA), length(pvalues[[1]])))
    
    # transform it to a DataFrame and return it
    return(.asDataFrame(pvalues_full, dataset@settings@sampleData[,sampleID]))
}


#'
#'
#' @noRd
.testPsiWithFisherPerGroup <- function(dataset, groupID, rawCounts, rawOtherCounts, internBPPARAM){
    # get group to test
    group <- sampleGroup(dataset@settings)
    group2Test <- group == groupID
    group2Test[is.na(group2Test)] <- FALSE
    
    fullFisherTable <- data.table(
        TP=rowSums(rawCounts[     , group2Test,with=FALSE]),
        FP=rowSums(rawCounts[     ,!group2Test,with=FALSE]),
        FN=rowSums(rawOtherCounts[, group2Test,with=FALSE]),
        TN=rowSums(rawOtherCounts[,!group2Test,with=FALSE])
    )
    
    # test only where at least the group has one read
    fisherTableToTest <- fullFisherTable[TP+FN > 0]
    pvalues <- unlist(bplapply(1:nrow(fisherTableToTest), BPPARAM = internBPPARAM, fisherTableToTest=fisherTableToTest,
            function(idx, fisherTableToTest){
                fisher.test(matrix(as.integer(fisherTableToTest[idx]), nrow=2))$p.value 
            }
    ))
    
    # add NAs wher the test group did not had any read
    fullFisherTable[,pvalue:=as.numeric(NA)]
    fullFisherTable[TP+FN>0,pvalue:=pvalues]
    return(fullFisherTable[,pvalue])
}



