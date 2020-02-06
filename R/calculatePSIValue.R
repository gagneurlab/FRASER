##
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## This file contains all functions for calculating the PSI values
## It calculates the PSI value for the junctions and
## the sitePSI value for intron retention
##

#'
#' PSI value calculation
#' 
#' This function calculates the PSI values for each junction and splice site
#' based on the FraseRDataSet object
#'
#' @inheritParams countRNA
#' @param types A vector with the psi types which should be calculated. Default 
#' is all of psi5, psi3 and psiSite.
#' @param overwriteCts FALSE or TRUE (the default) the total counts (aka N) will
#'              be recalculated based on the existing junction counts (aka K)
#' @return FraseRDataSet
#' @export
#' @examples
#'   fds <- makeExampleFraseRDataSet()
#'   fds <- calculatePSIValues(fds)
calculatePSIValues <- function(fds, types=psiTypes, overwriteCts=FALSE, 
                                BPPARAM=bpparam()){
    # check input
    stopifnot(is(fds, "FraseRDataSet"))
    
    # calculate PSI value for each sample
    for(psiType in unique(vapply(types, whichReadType, fds=fds, ""))){
        fds <- calculatePSIValuePrimeSite(fds, psiType=psiType,
                                overwriteCts=overwriteCts, BPPARAM=BPPARAM)
    }
    
    # calculate the delta psi value
    for(psiType in types){
        assayName <- paste0("delta_", psiType)
        fds <- calculateDeltaPsiValue(fds, psiType, assayName)
    }
    
    # save final Fraser object to disk
    fds <- saveFraseRDataSet(fds)
    
    # return it
    return(fds)
}


#'
#' calculates the PSI value for the given prime site of the junction
#'
#' @noRd
calculatePSIValuePrimeSite <- function(fds, psiType, overwriteCts, BPPARAM){
    stopifnot(is(fds, "FraseRDataSet"))
    stopifnot(isScalarCharacter(psiType))
    stopifnot(psiType %in% c("j", "ss"))
    
    if(psiType=="ss"){
        return(calculateSitePSIValue(fds, overwriteCts, BPPARAM=BPPARAM))
    }
    
    message(date(), ": Calculate the PSI 5 and 3 values ...")
    
    # generate a data.table from granges
    countData <- as.data.table(granges(rowRanges(fds, type=psiType)))
    
    # check if we have to compute N
    if(!all(paste0("rawOtherCounts_psi", c(5, 3)) %in% assayNames(fds))){
        overwriteCts <- TRUE
    }
    
    # calculate psi value
    psiValues <- bplapply(samples(fds), countData=countData,
            overwriteCts=overwriteCts, BPPARAM=BPPARAM,
        FUN=function(sample, countData, overwriteCts){
            
            # add sample specific counts to the data.table (K)
            sample <- as.character(sample)
            countData[,k:=list(K(fds, type="psi5")[,sample])]
            
            # get other counts (aka N) from cache or compute it
            if(isFALSE(overwriteCts)){
                countData[,o5:=counts(fds, type="psi5", side="oth")[,sample]]
                countData[,o3:=counts(fds, type="psi3", side="oth")[,sample]]
            } else {
                # compute other counts in strand specific way (+ and *) | (-)
                countData[,c("o5", "o3"):=list(0L, 0L)]
                plus <- countData[,strand %in% c("+", "*")]
                
                # compute psi5/3 on strand + and *
                countData[plus, o5:=sum(k)-k, by="seqnames,start"]
                countData[plus, o3:=sum(k)-k, by="seqnames,end"]
                
                # compute psi5/3 on strand -
                countData[!plus, o5:=sum(k)-k, by="seqnames,end"]
                countData[!plus, o3:=sum(k)-k, by="seqnames,start"]
            }
            
            # calculate psi value
            countData[,c("psi5", "psi3"):=list(k/(k+o5), k/(k+o3))]
            
            # if psi is NA this means there were no reads at all so set it to 1
            countData[is.na(psi5),psi5:=1]
            countData[is.na(psi3),psi3:=1]
            
            # write other counts and psi values to h5 file
            # get defind chunk sizes
            chunkDims <- c(
                min(nrow(countData), options()[['FraseR-hdf5-chunk-nrow']]),
                min(4, options()[['FraseR-hdf5-chunk-ncol']]))
            cacheFile <- getOtherCountsCacheFile(sample, fds)
            if(file.exists(cacheFile)){
                unlink(cacheFile)
            }
            writeHDF5Array(as.matrix(countData[,.(o5,o3,psi3,psi5)]), 
                            filepath=cacheFile, name="o5_o3_psi5_psi3", 
                            chunkdim=chunkDims, level=7, verbose=FALSE)
            
            # get counts as DelayedMatrix
            HDF5Array(filepath=cacheFile, name="o5_o3_psi5_psi3")
            
        }
    )
    names(psiValues) <- samples(fds)
    
    # merge it and assign it to our object
    assay(fds, type="j", "psi5") <- 
        saveAsHDF5(fds, "psi5", 
                    do.call(cbind, lapply(psiValues, function(x) x[,4])) )
    assay(fds, type="j", "psi3") <- 
        saveAsHDF5(fds, "psi3", 
                    do.call(cbind, lapply(psiValues, function(x) x[,3])) )
    
    if(isTRUE(overwriteCts)){
        assay(fds, type="j", "rawOtherCounts_psi5") <-
            saveAsHDF5(fds, "rawOtherCounts_psi5", 
                        do.call(cbind, lapply(psiValues, 
                                                function(x){ 
                                                    as.integer(x[,1])
                                                })) )
        assay(fds, type="j", "rawOtherCounts_psi3") <-
            saveAsHDF5(fds, "rawOtherCounts_psi3", 
                        do.call(cbind, lapply(psiValues, 
                                                function(x){
                                                    as.integer(x[,2])
                                                })) )
    }
    
    return(fds)
}


#'
#' This function calculates the site PSI values for each splice site
#' based on the FraseRDataSet object
#'
#' @noRd
calculateSitePSIValue <- function(fds, overwriteCts, BPPARAM){
    
    # check input
    stopifnot(is(fds, "FraseRDataSet"))
    
    message(date(), ": Calculate the PSI site values ...")
    
    psiName <- "psiSite"
    psiROCName <- "rawOtherCounts_psiSite"
    if(!psiROCName %in% assayNames(fds)){
        overwriteCts <- TRUE
    }
    
    # prepare data table for calculating the psi value
    countData <- data.table(
        spliceSiteID=c(
            rowData(fds, type="j")[["startID"]],
            rowData(fds, type="j")[["endID"]],
            rowData(fds, type="ss")[["spliceSiteID"]]
        ),
        type=rep(
            c("junction", "spliceSite"),
            c(length(fds)*2, length(nonSplicedReads(fds)))
        )
    )
    
    psiSiteValues <- bplapply(samples(fds), countData=countData, fds=fds,
        BPPARAM=BPPARAM, FUN=function(sample, countData, fds){
            if(verbose(fds) > 3){
                message("sample: ", sample)
            }
            # add sample specific counts to the data.table
            sample <- as.character(sample)
            sdata <- data.table(k=c(
                    rep(K(fds, type="psi3")[,sample], 2),
                    K(fds, type="psiSite")[,sample]))
            sdata <- cbind(countData, sdata)
            sdata[,os:=sum(k)-k, by="spliceSiteID"]
            
            # remove the junction part since we only want to calculate the
            # psi values for the splice sites themselfs
            sdata <- sdata[type=="spliceSite"]
            
            # calculate psi value
            sdata[,psiValue:=k/(os + k)]
            
            # if psi is NA this means there were no reads at all so set it to 1
            sdata[is.na(psiValue),psiValue:=1]
            
            # write other counts and psi values to h5 file
            # get defind chunk sizes
            chunkDims <- c(
                    min(nrow(sdata), options()[['FraseR-hdf5-chunk-nrow']]),
                    min(2, options()[['FraseR-hdf5-chunk-ncol']]))
            cacheFile <- getOtherCountsCacheFile(sample, fds)
            writeHDF5Array(as.matrix(sdata[,.(os, psiValue)]), 
                            filepath=cacheFile, name="oSite_psiSite", 
                            chunkdim=chunkDims, level=7, verbose=FALSE)
            
            # get counts as DelayedMatrix
            HDF5Array(filepath=cacheFile, name="oSite_psiSite")
        }
    )
    names(psiSiteValues) <- samples(fds)
    
    # merge it and assign it to our object
    assay(fds, type="ss", psiName) <- 
        saveAsHDF5(fds, psiName, 
                    do.call(cbind, lapply(psiSiteValues, function(x) x[,2])) )
    if(isTRUE(overwriteCts)){
        assay(fds, type="ss", psiROCName) <- 
            saveAsHDF5(fds, psiROCName, 
                    do.call(cbind, lapply(psiSiteValues, 
                                            function(x) {
                                                as.integer(x[,1])
                                            })) )
    }
    
    return(fds)
}

#'
#' calculates the delta psi value and stores it as an assay
#' @noRd
calculateDeltaPsiValue <- function(fds, psiType, assayName){
    
    message(date(), ": Calculate the delta for ", psiType, " values ...")
    
    # get psi values
    psiVal <- as.matrix(assays(fds)[[psiType]])
    
    # psi - median(psi)
    rowmedian <- rowMedians(psiVal, na.rm = TRUE)
    deltaPsi  <- psiVal - rowmedian
    
    # use as.matrix to rewrite it as a new hdf5 array
    assays(fds, type=psiType)[[assayName]] <- deltaPsi
    
    return(fds)
}

#'
#' returns the name of the cache file for the given sample
#' @noRd
getOtherCountsCacheFile <- function(sampleID, fds){
    # cache folder
    cachedir <- getOtherCountsCacheFolder(fds)
    
    # file name
    filename <- paste0("otherCounts-", sampleID, ".h5")
    
    # return it
    return(file.path(cachedir, filename))
}

#'
#' returns the name of the cache folder if caching is enabled 
#' @noRd
getOtherCountsCacheFolder <- function(fds){
    
    # cache folder
    cachedir <- file.path(workingDir(fds), "cache", "otherCounts", 
                            nameNoSpace(name(fds)))
    if(!dir.exists(cachedir)){
        dir.create(cachedir, recursive=TRUE)
    }
    
    # return it
    return(cachedir)
}