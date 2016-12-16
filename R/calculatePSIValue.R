##
## @author Christian Mertes \email{mertes@@in.tum.de}
## 
## This file contains all functions for calculating the PSI values
## It calculates the PSI value for the junctions and
## the sitePSI value for intron retention
##

#'
#' This function calculates the PSI values for each junction
#' based on the FraseRDataSet object
#' 
#' @export
#' 
calculatePSIValues <- function(dataset){
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")
    
    # generate a data.table from granges
    rawCountData <- dataset@splitReads
    countData <- cbind(
        data.table(
            chr = as.factor(seqnames(rawCountData)),
            start = start(rawCountData),
            end = end(rawCountData),
            strand = as.factor(strand(rawCountData)),
            counts = NA
        ),
        as.data.table(assays(rawCountData)$rawCounts)
    )
    
    # calculate 3/5' PSI for each sample
    assays(rawCountData)$psi3 <- .calculatePSIValuePrimeSite(
            countData=countData, psiType="3'",
            settings=dataset@settings
    )
    assays(rawCountData)$psi5 <- .calculatePSIValuePrimeSite(
            countData=countData, psiType="5'",
            settings=dataset@settings
    )
    
    dataset@splitReads <- rawCountData
    return(dataset)
}


#'
#' calculates the PSI value for the given prime site of the junction
#' 
.calculatePSIValuePrimeSite <- function(countData, settings, psiType){
    
    # convert psi type to the position of interest 
    psiCol <- ifelse(psiType == "3'", "start", "end")
    
    # calculate psi value
    psiValues <- bplapply(settings@sampleData[,sampleID],
                          countData=countData, psiCol=psiCol,
                          BPPARAM=settings@parallel, 
                           FUN = function(sample, countData, psiCol){
                               suppressPackageStartupMessages(library(FraseR))
                               
                               # check name, due to conversion
                               if(grepl("^\\d+$", sample)){
                                   sample <- paste0("X", sample)
                               }
                               
                               # init psi
                               countData[,psiValue:=as.numeric(NA)]
                               
                               # calculate psi value
                               countData[,psiValue:=get(sample)/sum(get(sample), na.rm = TRUE),
                                    by=eval(paste0("chr,", psiCol, ",strand"))
                               ]
                               
                               return(countData[,psiValue])
                           }
    )
    
    # merge it and set the column names
    df = DataFrame(matrix(unlist(psiValues), ncol = nrow(settings@sampleData)))
    names(df) <- settings@sampleData[,sampleID] 
    
    return(df)
}


#'
#' This function calculates the site PSI values for each splice site
#' based on the FraseRDataSet object
#' 
#' @export
#' 
calculateSitePSIValue <- function(dataset){
    
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")
    
    # get only the junctions
    splitReads <- dataset@splitReads
    nonSplicedCounts <- dataset@nonSplicedReads
    
    
    for(spliceType in c("Acceptor", "Donor")){
        modified_rows <- rowData(nonSplicedCounts)$type == spliceType
        splice_site   <- nonSplicedCounts[modified_rows]
        splice_ranges <- rowRanges(splice_site)
        
        # shift for start/end overlap
        splice_ranges <- GenomicRanges::shift(splice_ranges, ifelse(spliceType == "Acceptor", -1, 1))
        
        # find overlap
        overlap <- findOverlaps(splice_ranges, 
                splitReads, type=ifelse(spliceType == "Acceptor", "end", "start")
        )
        
        # sum up the junctions per site per sample
        junction_counts_per_site <- bplapply(
                dataset@settings@sampleData[,sampleID], 
                overlap=overlap, splitReads=splitReads,
                BPPARAM=dataset@settings@parallel,
                FUN = function(sample, overlap, splitReads){
                        suppressPackageStartupMessages(library(FraseR))
                    
                        dt <- data.table(
                                from = overlap@from,
                                count = as.numeric(assays(splitReads[,sample])$rawCounts[,1])[overlap@to]
                        )
                        dt <- dt[,.(counts=sum(count, na.rm=TRUE), is.na=all(is.na(count))), by=from]
                        return(dt[,counts])
                }
        )
        
        # add junction counts
        splice_site <- .mergeAssay("rawSplitReadCounts", splice_site,
                junction_counts_per_site, unique(overlap@from)
        )
       
        # wrote it to the data
        dataset@nonSplicedReads <- .mergeAssay("rawSplitReadCounts", dataset@nonSplicedReads,
                assays(splice_site)$rawSplitReadCounts, modified_rows
        )
    }
    
    
    # add junction persent spliced in value (sitePSI) 
    sitePSIValue <- .getAssayAsDataTable(dataset@nonSplicedReads, "rawSplitReadCounts") / (
        .getAssayAsDataTable(dataset@nonSplicedReads, "rawSplitReadCounts") +
            .getAssayAsDataTable(dataset@nonSplicedReads, "rawCounts")
    )
    dataset@nonSplicedReads <- .mergeAssay("sitePSI", dataset@nonSplicedReads, sitePSIValue)
    
    return(dataset)
}


#'
#' convert a data.table to a DataFrame and keep the colnames
#' 
.asDataFrame <- function(dataframe, colname = colnames(dataframe)){
    dataframe <- DataFrame(dataframe)
    colnames(dataframe) <- colname
    return(dataframe)
}


#'
#' @param assay name of assay to use
#' @param se summarizedExperiment object to add the assay
#' @param df dataframe/data to add as assay
#' 
.mergeAssay <- function(assay, se1, df1, rowidx1 = 1:dim(se1)[1]){
    # create empty dataframe for assay if not present yet
    if(!any(names(assays(se1)) %in% assay)){
        tmp_df <- DataFrame(matrix(as.numeric(NA), nrow = dim(se1)[1], ncol = dim(se1)[2]))
        colnames(tmp_df) <- colnames(se1)
        assays(se1)[[assay]] <- tmp_df
    } 
    
    # make a data frame out of the given data if needed
    if(!any(class(df1) %in% "DataFrame")){
        df1 <- .asDataFrame(df1, colnames(se1))
    }
    
    # merge the data
    tmp_df <- assays(se1)[[assay]]
    tmp_df[rowidx1,] <- df1
    assays(se1)[[assay]] <- tmp_df
    
    # return it
    return(se1)
}


#'
#' get the assay as data.table from a SummarizedExperiment object
#' 
.getAssayAsDataTable <- function(se, assay, na_as_zero = TRUE){
    if(!any(names(assays(se)) %in% assay)){
        stop("The given assay: '", assay, "' is not present in this object")
    }
    dt <- as.data.table(assays(se)[[assay]])
    if(na_as_zero){
        dt[is.na(dt)] <- 0
    }
    colnames(dt) <- colnames(assays(se)[[assay]])
    return(dt)
}

