##
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for calculating the statistics
## First it starts with calculating the Z-score for each
## site and then the p-values are calculated dependend on the
## given method in the setting file
##

#'
#' Calculate the zscore for each PSI value.
#'
#' @export
#' @examples
#'   fds <- countRNAData(createTestFraseRSettings())
#'   fds <- calculatePSIValues(fds)
#'   fds <- calculateZScores(fds)
calculateZScores <- function(fds, type=psiTypes){

    # check input
    stopifnot(class(fds) == "FraseRDataSet")

    # calculate zscore for each psi type
    for(pt in type){
        fds <- calculateZScorePerDataSet(fds, pt)
    }

    return(fds)
}

#'
#' calculates the zscore for a given data type and a given psi type
#' and adds it directly to the dataset itself
#'
#' @noRd
calculateZScorePerDataSet <- function(fds, psiType){

    message(date(), ": Calculate the Zscore for ", psiType, " values ...")

    # get raw data and replace NA's with zeros
    psiVal <- assay(fds, psiType)

    # z = ( x - mean ) / sd
    rowmean <- rowMeans(psiVal, na.rm = TRUE)
    rowsd   <- apply(psiVal, 1, sd, na.rm = TRUE)
    zscores <- (psiVal - rowmean) / rowsd

    # use as.matrix to rewrite it as a new hdf5 array
    zScores(fds, type=psiType) <- as.matrix(zscores)

    return(fds)
}

#'
#' calculates the P-Value for the given FraseR dataset object
#' The P-Value calculation is based on the given method in the
#' FraseRSettings object
#'
#' @export
#' @examples
#'   fds <- countRNAData(createTestFraseRSettings())
#'   fds <- calculatePSIValues(fds)
#'   fds <- calculatePValues(fds)
calculatePValues <- function(fds, type=psiTypes, internBPPARAM=bpparam(), ...){
    # check input
    stopifnot(class(fds) == "FraseRDataSet")
    enforceHDF5 <- FALSE

    # get correct method:
    FUN <- switch(method(fds),
        betaBin = {
            enforceHDF5 <- TRUE
            list(FUN=pvalueByBetaBinomialPerType, pvalFun=betabinVglmTest)
        },
        Fisher  = {
            list(FUN=pvalueByFisherPerType, internBPPARAM=internBPPARAM)
        },
        stop("The provided method is not present for this package.",
             "Please set the method to one of the following:",
             "Fisher, betaBin, DESeq2, Martin"
        )
    )


    # check, that the object is stored as HDF5 array if requested
    aIsHDF5 <- sapply(assays(fds), function(x) any("DelayedArray" == is(x)))
    if(enforceHDF5 & !all(aIsHDF5)){
        message(date(), ": The data is not stored in a HDF5Array. ",
                "To improve the performance we will store now ",
                "the data in HDF5 format.")
        # fds <- saveFraseRDataSet(fds)
    }

    # test all 3 different types
    for(psiType in type){
        fds <- do.call(FUN[[1]],
            c(fds=fds, aname=aname, psiType=psiType, FUN[-1], ...)
        )
        fds <- saveFraseRDataSet(fds)
        gc()
    }

    # return the new datasets
    return(fds)
}


#'
#' calculates the pvalue per type (psi3,psi5,spliceSite) with fisher
#'
#' @noRd
pvalueByFisherPerType <- function(dataset, psiType, internBPPARAM){
    # go over each group but no NA's
    group   <- sampleGroup(dataset)

    # reads to test for abberent splicing (eg: nonSplicedReads)
    rawCounts <- counts(dataset, type=psiType, side="ofIn")
    rawOtherCounts <- counts(dataset, type=psiType, side="other")

    pvalues <- bplapply(unique(na.omit(group)), dataset=dataset,
            rawCounts=rawCounts, rawOtherCounts=rawOtherCounts,
            BPPARAM=parallel(dataset),
            internBPPARAM=internBPPARAM,
            FUN=.testPsiWithFisherPerGroup
    )
    names(pvalues) <- as.character(unique(na.omit(group)))
    pvalues_full <- pvalues[as.character(group)]

    # add NA's to the non tested ones
    pvalues_full[is.na(group)] <- list(
            rep(as.numeric(NA), length(pvalues[[1]]))
    )

    # transform it to a DataFrame and return it
    return(.asDataFrame(pvalues_full, samples(dataset)))
}


#'
#' calculates the pvalue per group with fisher
#'
#' @noRd
.testPsiWithFisherPerGroup <- function(dataset, groupID, rawCounts,
            rawOtherCounts, internBPPARAM){
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
    pvalues <- unlist(bplapply(1:nrow(fisherTableToTest),
            BPPARAM=internBPPARAM,
            fisherTableToTest=fisherTableToTest,
            FUN=function(idx, fisherTableToTest){
                fisher.test(matrix(as.integer(
                    fisherTableToTest[idx]), nrow=2
                ))$p.value
            }
    ))

    # add NAs wher the test group did not had any read
    fullFisherTable[,pvalue:=as.numeric(NA)]
    fullFisherTable[TP+FN>0,pvalue:=pvalues]
    return(fullFisherTable[,pvalue])
}


