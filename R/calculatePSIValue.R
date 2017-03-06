##
## @author Christian Mertes \email{mertes@@in.tum.de}
## 
## This file contains all functions for calculating the PSI values
## It calculates the PSI value for the junctions and
## the sitePSI value for intron retention
##

#'
#' This function calculates the PSI values for each junction and splice site
#' based on the FraseRDataSet object
#' 
#' @export
#' @examples
#'   fds <- counRNAData(createTestFraseRSettings())
#'   fds <- calculatePSIValues(fds)
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
    message(date(), ": Calculate the PSI3 values ...")
    psiValues <- .calculatePSIValuePrimeSite(
            countData=countData, psiType="3'",
            settings=dataset@settings
    )
    assays(rawCountData)$psi3 <- psiValues[["psi"]]
    assays(rawCountData)$rawOtherCounts_psi3 <- psiValues[["counts"]]
    
    message(date(), ": Calculate the PSI5 values ...")
    psiValues <- .calculatePSIValuePrimeSite(
            countData=countData, psiType="5'",
            settings=dataset@settings
    )
    assays(rawCountData)$psi5 <- psiValues[["psi"]]
    assays(rawCountData)$rawOtherCounts_psi5 <- psiValues[["counts"]]
    
    dataset@splitReads <- rawCountData
    
    # calculate siteSplice values
    dataset <- .calculateSitePSIValue(dataset)
    
    # return it
    return(dataset)
}


#'
#' calculates the PSI value for the given prime site of the junction
#' 
#' @noRd
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
                           
            # calculate other split read counts
            countData[,rawOtherCounts:=sum(get(sample), na.rm=TRUE) - get(sample), 
                    by=eval(paste0("chr,", psiCol, ",strand"))
            ]
                           
            # calculate psi value
            countData[,psiValue:=get(sample)/(get(sample) + rawOtherCounts)]
                           
            return(list(
                    rawOtherCounts=countData[,rawOtherCounts],
                    psiValue=countData[,psiValue]
            ))
        }
    )
    
    # merge it and set the column names
    dfPsi <- DataFrame(matrix(unlist(sapply(psiValues, "[", "psiValue")), 
            ncol = nrow(settings@sampleData)
    ))
    names(dfPsi) <- settings@sampleData[,sampleID] 
    dfCounts  <- DataFrame(matrix(unlist(sapply(psiValues, "[", "rawOtherCounts")), 
            ncol = nrow(settings@sampleData)
    ))
    names(dfCounts) <- settings@sampleData[,sampleID] 
    
    return(SimpleList(psi=dfPsi, counts=dfCounts))
}


#'
#' This function calculates the site PSI values for each splice site
#' based on the FraseRDataSet object
#' 
#' @noRd
.calculateSitePSIValue <- function(dataset){
    
    # check input
    stopifnot(class(dataset) == "FraseRDataSet")
    
    # get only the junctions
    splitReads <- dataset@splitReads
    nonSplicedCounts <- dataset@nonSplicedReads
    
    
    for(spliceType in c("Acceptor", "Donor")){
        message(date(), ": Calculate the splice site values for ", spliceType, " ...")
        modified_rows <- rowData(nonSplicedCounts)$type == spliceType
        splice_site   <- nonSplicedCounts[modified_rows]
        splice_ranges <- rowRanges(splice_site)
        
        # shift for start/end overlap
        splice_ranges <- shift(splice_ranges, 
                ifelse(spliceType == "Acceptor", -1, 1)
        )
        
        # find overlap
        overlap <- findOverlaps(splice_ranges, splitReads, 
                type=ifelse(spliceType == "Acceptor", "end", "start")
        )
        
        # sum up the junctions per site per sample
        junction_counts_per_site <- bplapply(
                    samples(dataset@settings), 
                    overlap=overlap, splitReads=splitReads,
                    BPPARAM=parallel(dataset@settings),
                    FUN = function(sample, overlap, splitReads){
            suppressPackageStartupMessages(library(FraseR))
        
            dt <- data.table(
                from = overlap@from,
                count = as.numeric(assays(
                        splitReads[,sample])$rawCounts[,1])[overlap@to]
            )
            dt <- dt[,.(counts=sum(count, na.rm=TRUE), 
                    is.na=all(is.na(count))), by=from
            ]
            return(dt[,counts])
        })
        
        # add junction counts
        splice_site <- .mergeAssay("rawSplitReadCounts", splice_site,
                junction_counts_per_site, unique(overlap@from)
        )
       
        # wrote it to the data
        dataset@nonSplicedReads <- .mergeAssay("rawOtherCounts_sitePSI", 
                dataset@nonSplicedReads,
                assays(splice_site)$rawSplitReadCounts, modified_rows
        )
    }
    
    
    # add junction persent spliced in value (sitePSI) 
    sitePSIValue <- 
        .getAssayAsDataTable(dataset@nonSplicedReads, "rawOtherCounts_sitePSI") / (
            .getAssayAsDataTable(dataset@nonSplicedReads, "rawOtherCounts_sitePSI") +
            .getAssayAsDataTable(dataset@nonSplicedReads, "rawCounts")
    )
    dataset@nonSplicedReads <- FraseR:::.mergeAssay("sitePSI", dataset@nonSplicedReads, sitePSIValue)
    
    return(dataset)
}


#'
#' merge assays into one
#' 
#' @param assay name of assay to use
#' @param se summarizedExperiment object to add the assay
#' @param df dataframe/data to add as assay
#' 
#' @noRd
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



