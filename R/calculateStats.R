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
    
    dataset <- FraseR:::.calculateZScorePerDataSet(dataset, "splitReads", "psi3")
    dataset <- FraseR:::.calculateZScorePerDataSet(dataset, "splitReads", "psi5")
    dataset <- FraseR:::.calculateZScorePerDataSet(dataset, "nonSplicedReads", "sitePSI")
    
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
    psiVal <- FraseR:::.getAssayAsDataTable(seCounts, psiType)
    
    # z = ( x - mean ) / sd
    rowmean <- rowMeans(psiVal, na.rm = TRUE)
    rowsd   <- apply(psiVal, 1, sd, na.rm = TRUE)
    zscores <- (psiVal - rowmean) / rowsd
    
    # add it to the FraseR object
    assayName <- paste0("zscore_", psiType)
    zscores   <- FraseR:::.asDataFrame(zscores, dataset@settings@sampleData[,sampleID])
    assays(slot(dataset, readType))[[assayName]] <- zscores
          
    return(dataset)  
}

#'
#' calculates the P-Value for the given FraseR dataset object
#' The P-Value calculation is based on the given method in the
#' FraseRSettings object
#'
#' @export
calculatePValues <- function(dataset){
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")
    
    # check which method we should use
    method <- dataset@settings@method
    if(method == "Fisher"){
        return(FraseR:::.testPsiWithFisher(dataset))
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
.testPsiWithFisher <- function(dataset){
    
    # test all 3 different types
    assays(dataset@splitReads)$pvalue_psi3 <- 
            .testPsiWithFisherPerType(dataset, "splitReads", "psi3")
    assays(dataset@splitReads)$pvalue_psi5 <- 
            .testPsiWithFisherPerType(dataset, "splitReads", "psi5")
    assays(dataset@nonSplicedReads)$pvalue_sitePSI <- 
            .testPsiWithFisherPerType(dataset, "nonSplicedReads", "sitePSI")
    
    # return the new datasets
    return(dataset)
}

#'
#'
#' @noRd
.testPsiWithFisherPerType <- function(dataset, readType, psiType){
    # go over each group but no NA's
    group   <- sampleGroup(dataset@settings)
    
    # reads to test for abberent splicing (eg: nonSplicedReads)
    rawCounts <- FraseR:::.getAssayAsDataTable(slot(dataset, readType), "rawCounts")
    
    # other reads (eg: splitReads)
    rawOtherCounts <- FraseR:::.getAssayAsDataTable(slot(dataset, readType), paste0("rawOtherCounts_", psiType))
    
    pvalues <- bplapply(unique(na.omit(group)), dataset=dataset, 
                        rawCounts=rawCounts, rawOtherCounts=rawOtherCounts,
                        BPPARAM=dataset@settings@parallel,
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
.testPsiWithFisherPerGroup <- function(dataset, groupID, rawCounts, rawOtherCounts){
    # get group to test
    group <- sampleGroup(dataset@settings)
    group2Test <- group == groupID
    group2Test[is.na(group2Test)] <- FALSE
    
    # groups to test against 
    pvalues <- sapply(1:dim(rawCounts)[1], function(idx){
        fisher.test(matrix(nrow=2,
                           c(
                               TP=sum(rawCounts[     idx, group2Test,with=FALSE]),
                               FP=sum(rawCounts[     idx,!group2Test,with=FALSE]),
                               FN=sum(rawOtherCounts[idx, group2Test,with=FALSE]),
                               TN=sum(rawOtherCounts[idx,!group2Test,with=FALSE])
                           )
        ))$p.value
    })
    return(pvalues)
}


#' 
convert_dataframe_columns_to_Rle <- function(data, index2convert = 1:dim(dataframe)[2]){
    
    # convert all given indices
    for (i in index2convert){
        data[,i] <- Rle(data[,i])
    }
    
    return(data)
}




#'
#' Filter the data based on a minimum of expression level over all samples
#' It removes the junction and also the corresponding Donor and Acceptor site
#' within the SummarizedExperiment object
#' @noRd
filter_junction_data <- function(data, minExpRatio = 0.8){
    # get only the junctions
    junctions <- which(rowData(data)$type == "Junction")
    
    # get the expression counts for each junction
    dt <- get_assay_as_data_table(data, "counts", FALSE)[junctions]
    
    # calculate the expression ratio per site
    expression <- apply(dt, 1, function(x) sum(!(is.na(x) | x == 0)))
    expression <- expression / dim(dt)[2]
    
    cutoff <- expression >= minExpRatio
    
    # get the hits (junction/acceptor/donor) in our full data set
    hits <- unique(unlist(sapply(c("start", "end"), function(type){
        findOverlaps(type = type,
                     rowRanges(data)[junctions][cutoff], 
                     shift(rowRanges(data), ifelse(type == "start", 1, -1))
        )@to
    })))
    junction_sites <- hits[rowData(data)$type[hits] != "Junction"]
    
    # filter the object and return it
    return(data[c(junction_sites, junctions[cutoff])])
}
