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
    
    dataset <- .calculateZScorePerDataSet(dataset, "splitReads", "psi3")
    dataset <- .calculateZScorePerDataSet(dataset, "splitReads", "psi5")
    dataset <- .calculateZScorePerDataSet(dataset, "nonSplicedReads", "sitePSI")
    
    return(dataset)
}

#'
#' calculates the zscore for a given data type and a given psi type
#' and adds it directly to the dataset itself
#'
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
#' 
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
