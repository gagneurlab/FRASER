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
#' based on the FraserDataSet object
#'
#' @inheritParams countRNA
#' @param types A vector with the psi types which should be calculated. Default 
#' is all of psi5, psi3 and theta.
#' @param overwriteCts FALSE or TRUE (the default) the total counts (aka N) will
#'              be recalculated based on the existing junction counts (aka K)
#' @return FraserDataSet
#' @export
#' @examples
#'   fds <- createTestFraserDataSet()
#'   fds <- calculatePSIValues(fds, types="psi5")
#'   
#'   ### usually one would run this function for all psi types by using:
#'   # fds <- calculatePSIValues(fds)
calculatePSIValues <- function(fds, types=psiTypes, overwriteCts=FALSE, 
                                BPPARAM=bpparam()){
    # check input
    stopifnot(is(fds, "FraserDataSet"))
    
    # calculate PSI value for each sample
    for(psiType in unique(vapply(types, whichReadType, fds=fds, ""))){
        fds <- calculatePSIValuePrimeSite(fds, psiType=psiType,
                                overwriteCts=overwriteCts, BPPARAM=BPPARAM)
    }
    
    # calculate intron jaccard index
    fds <- calculateIntronJaccardIndex(fds)
    
    # calculate the delta psi value
    for(psiType in types){
        assayName <- paste0("delta_", psiType)
        fds <- calculateDeltaPsiValue(fds, psiType, assayName)
    }
    
    # return it
    return(fds)
}


#'
#' calculates the PSI value for the given prime site of the junction
#'
#' @noRd
calculatePSIValuePrimeSite <- function(fds, psiType, overwriteCts, BPPARAM){
    stopifnot(is(fds, "FraserDataSet"))
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
    
    h5DatasetName <- "o5_o3_psi5_psi3"
    
    # calculate psi value
    psiValues <- bplapply(samples(fds), countData=countData,
            overwriteCts=overwriteCts, BPPARAM=BPPARAM,
        FUN=function(sample, countData, overwriteCts){
            
            # get sample
            sample <- as.character(sample)
            
            # check if other counts and psi values chache file exists already
            cacheFile <- getOtherCountsCacheFile(sample, fds)
            ans <- checkPsiCacheFile(cFile=cacheFile, dName=h5DatasetName, 
                    overwrite=overwriteCts, ptype=c("psi5", "psi3"), fds=fds)
            if(!is.null(ans)){
                return(ans)
            }
            
            # add sample specific counts to the data.table (K)
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
            
            # if no HDF5 is requested return it as matrix
            if(dontWriteHDF5(fds)){
                return(DelayedArray(as.matrix(countData[,.(o5,o3,psi5,psi3)])))
            }
            
            # write other counts and psi values to h5 file
            # get defined chunk sizes
            chunkDims <- c(
                min(nrow(countData), options()[['FRASER-hdf5-chunk-nrow']]),
                1)
            writeHDF5Array(as.matrix(countData[,.(o5,o3,psi5,psi3)]), 
                            filepath=cacheFile, name=h5DatasetName,
                            chunkdim=chunkDims, level=7, verbose=FALSE)
            
            # get counts as DelayedMatrix
            HDF5Array(filepath=cacheFile, name=h5DatasetName)
        }
    )
    names(psiValues) <- samples(fds)
    
    # merge it and assign it to our object
    assay(fds, type="j", "psi5", withDimnames=FALSE) <- do.call(cbind, 
            mapply('[', psiValues, TRUE, 3, drop=FALSE, SIMPLIFY=FALSE))
    assay(fds, type="j", "psi3", withDimnames=FALSE) <- do.call(cbind, 
            mapply('[', psiValues, TRUE, 4, drop=FALSE, SIMPLIFY=FALSE))
    
    if(isTRUE(overwriteCts)){
        assay(fds, type="j", "rawOtherCounts_psi5", withDimnames=FALSE) <- 
            do.call(cbind, mapply('[', psiValues, TRUE, 1,
                    drop=FALSE, SIMPLIFY=FALSE))
        assay(fds, type="j", "rawOtherCounts_psi3", withDimnames=FALSE) <- 
            do.call(cbind, mapply('[', psiValues, TRUE, 2,
                    drop=FALSE, SIMPLIFY=FALSE))
    }
    
    return(fds)
}


#'
#' Helper function to check for PSI value cached data
#'
#' @noRd
checkPsiCacheFile <- function(cFile, dName, overwrite, ptype, fds){
    if(file.exists(cFile) && dName %in% h5ls(cFile)$name){
        h5 <- HDF5Array(filepath=cFile, name=dName)
        aNames <- paste0("rawOtherCounts_", ptype)
        if(isFALSE(overwrite) && all( aNames %in% assayNames(fds)) &&
                nrow(h5) == nrow(K(fds, type=ptype[1]))){
            return(h5)
        }
        h5delete(cFile, name=dName)
    }
    return(NULL)
}


#'
#' This function calculates the site PSI values for each splice site
#' based on the FraserDataSet object
#'
#' @noRd
calculateSitePSIValue <- function(fds, overwriteCts, BPPARAM){
    
    # check input
    stopifnot(is(fds, "FraserDataSet"))
    
    message(date(), ": Calculate the PSI site values ...")
    
    psiName <- "theta"
    psiROCName <- "rawOtherCounts_theta"
    if(!psiROCName %in% assayNames(fds)){
        overwriteCts <- TRUE
    }
    psiH5datasetName <- "oSite_theta"
    
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
    
    thetaValues <- bplapply(samples(fds), countData=countData, fds=fds,
        BPPARAM=BPPARAM, FUN=function(sample, countData, fds){
            if(verbose(fds) > 3){
                message("sample: ", sample)
            }
            
            # get sample
            sample <- as.character(sample)
            
            # get counts and theta values from cache file if it exists
            cacheFile <- getOtherCountsCacheFile(sample, fds)
            ans <- checkPsiCacheFile(cFile=cacheFile, dName=psiH5datasetName, 
                    overwrite=overwriteCts, ptype="theta", fds=fds)
            if(!is.null(ans)){
                return(ans)
            }
            
            # add sample specific counts to the data.table
            sdata <- data.table(k=c(
                    rep(K(fds, type="psi3")[,sample], 2),
                    K(fds, type="theta")[,sample]))
            sdata <- cbind(countData, sdata)
            sdata[,os:=sum(k)-k, by="spliceSiteID"]
            
            # remove the junction part since we only want to calculate the
            # psi values for the splice sites themselves
            sdata <- sdata[type=="spliceSite"]
            
            # calculate psi value
            sdata[,psiValue:=k/(os + k)]
            
            # if psi is NA this means there were no reads at all so set it to 1
            sdata[is.na(psiValue),psiValue:=1]
            
            # if no HDF5 is requested return it as matrix
            if(dontWriteHDF5(fds)){
                return(DelayedArray(as.matrix(sdata[,.(os, psiValue)])))
            }
            
            # write other counts and psi values to h5 file
            # get defined chunk sizes
            chunkDims <- c(
                    min(nrow(sdata), options()[['FRASER-hdf5-chunk-nrow']]),
                    2)
            writeHDF5Array(as.matrix(sdata[,.(os, psiValue)]), 
                            filepath=cacheFile, name=psiH5datasetName, 
                            chunkdim=chunkDims, level=7, verbose=FALSE)
            
            # get counts as DelayedMatrix
            HDF5Array(filepath=cacheFile, name=psiH5datasetName)
        }
    )
    names(thetaValues) <- samples(fds)
    
    # merge it and assign it to our object
    assay(fds, type="ss", psiName, withDimnames=FALSE) <- do.call(cbind, 
            mapply('[', thetaValues, TRUE, 2, drop=FALSE, 
                    SIMPLIFY=FALSE))
    if(isTRUE(overwriteCts)){
        assay(fds, type="ss", psiROCName, withDimnames=FALSE) <- do.call(cbind,
                mapply('[', thetaValues, TRUE, 1, drop=FALSE, 
                        SIMPLIFY=FALSE))
    }
    
    return(fds)
}

#'
#' calculates the delta psi value and stores it as an assay
#' @noRd
calculateDeltaPsiValue <- function(fds, psiType, assayName){
    
    message(date(), ": Calculate the delta for ", psiType, " values ...")
    
    # get psi values
    psiVal <- assays(fds)[[psiType]]
    
    # psi - median(psi)
    rowmedian <- rowMedians(psiVal, na.rm = TRUE)
    deltaPsi  <- psiVal - rowmedian
    
    # rewrite it as a new hdf5 array
    assay(fds, assayName, type=psiType, withDimnames=FALSE) <- deltaPsi

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

#'
#' calculates the intron jaccard index values for all junctions
#' @inheritParams countRNA
#' 
calculateIntronJaccardIndex <- function(fds){
    # retrieve junction and splice site annotation with count information
    junction_dt <- data.table(
        as.data.table(rowRanges(fds, type="j"))[,.(seqnames, start, end, 
                                                   strand, startID, endID)], 
        as.matrix(K(fds, type="psi5")))
    ss_dt <- data.table(
        as.data.table(rowRanges(fds, type="ss"))[,.(seqnames, start, end, 
                                                    strand, spliceSiteID)], 
        as.matrix(N(fds, type="theta")))
    
    # melt to have one row per sample - junction combination
    junction_dt <- melt(junction_dt, variable.name="sampleID", value.name="k", 
                        id.vars=c("seqnames", "start", "end", "strand", 
                                  "startID", "endID"))
    ss_dt <- melt(ss_dt, variable.name="sampleID", value.name="n", 
                  id.vars=c("seqnames", "start", "end", "strand", 
                            "spliceSiteID"))
    
    # merge junction information with splice site annotation (theta)
    junction_dt <- merge(junction_dt, ss_dt[,.(spliceSiteID, sampleID, n)], 
                         all.x=TRUE, by.x=c("startID", "sampleID"), 
                         by.y=c("spliceSiteID", "sampleID"), sort=FALSE)
    setnames(junction_dt, "n", "n_donor")
    junction_dt <- merge(junction_dt, ss_dt[,.(spliceSiteID, sampleID, n)], 
                         all.x=TRUE, by.x=c("endID", "sampleID"), 
                         by.y=c("spliceSiteID", "sampleID"), sort=FALSE)
    setnames(junction_dt, "n", "n_acceptor")
    rm(ss_dt)
    gc()
    
    # TODO deal with missing endIDs (why do we have them? lost in filtering?)
    # for now: replace n_donor with N from psi5 (and same for n_acceptor and psi3)
    junction_dt_nacceptor <- data.table(
        as.data.table(rowRanges(fds, type="j"))[,.(startID, endID)], 
        as.matrix(N(fds, type="psi3")))
    junction_dt_nacceptor <- melt(junction_dt_nacceptor, 
                        variable.name="sampleID", value.name="n_psi3", 
                        id.vars=c("startID", "endID"))
    junction_dt[, n_psi3:=junction_dt_nacceptor[,n_psi3]]
    rm(junction_dt_nacceptor)
    gc()
    junction_dt_ndonor <- data.table(
        as.data.table(rowRanges(fds, type="j"))[,.(startID, endID)], 
        as.matrix(N(fds, type="psi5")))
    junction_dt_ndonor <- melt(junction_dt_ndonor, 
                               variable.name="sampleID", value.name="n_psi5", 
                               id.vars=c("startID", "endID"))
    junction_dt[, n_psi5:=junction_dt_ndonor[,n_psi5]]
    rm(junction_dt_ndonor)
    gc()
    
    # replace n (with non-split counts) by n_psi3/n_psi5 if NA
    junction_dt[is.na(n_acceptor), n_acceptor:=n_psi3]
    junction_dt[is.na(n_donor), n_donor:=n_psi5]
    junction_dt[n_acceptor < n_psi3, n_acceptor:=n_psi3]
    junction_dt[n_donor < n_psi5, n_donor:=n_psi5]
    
    # calculate intron_jaccard
    junction_dt[, denominator:=(n_donor + n_acceptor - k)]
    junction_dt[, intron_jaccard:= k / denominator]
    junction_dt[is.nan(intron_jaccard), intron_jaccard:=1] # n_donor = n_acceptor = 0
    
    # convert to matrix to store it as assay in the fds
    intron_jaccard <- matrix(junction_dt[,intron_jaccard], nrow=nrow(fds), ncol=ncol(fds), byrow=FALSE)
    rownames(intron_jaccard) <- rownames(fds)
    colnames(intron_jaccard) <- colnames(fds)
    assay(fds, type="j", "intron_jaccard", withDimnames=FALSE) <- intron_jaccard
    
    # store denominator
    jaccard_denom <- matrix(junction_dt[,denominator], nrow=nrow(fds), ncol=ncol(fds), byrow=FALSE)
    rownames(jaccard_denom) <- rownames(fds)
    colnames(jaccard_denom) <- colnames(fds)
    assay(fds, type="j", "rawOtherCounts_intron_jaccard", 
            withDimnames=FALSE) <- jaccard_denom - matrix(junction_dt[,k], 
                                                            nrow=nrow(fds), 
                                                            ncol=ncol(fds), 
                                                            byrow=FALSE)
    return(fds)
}

calculatePhi <- function(fds){
    # retrieve junction and splice site annotation with count information
    junction_dt <- data.table(
        as.data.table(rowRanges(fds, type="j"))[,.(seqnames, start, end, 
                                                   strand, startID, endID)], 
        as.matrix(K(fds, type="psi5")))
    ss_dt <- data.table(
        as.data.table(rowRanges(fds, type="ss"))[,.(seqnames, start, end, 
                                                    strand, spliceSiteID)], 
        as.matrix(N(fds, type="theta")))
    
    # melt to have one row per sample - junction combination
    junction_dt <- melt(junction_dt, variable.name="sampleID", value.name="k", 
                        id.vars=c("seqnames", "start", "end", "strand", 
                                  "startID", "endID"))
    ss_dt <- melt(ss_dt, variable.name="sampleID", value.name="n", 
                  id.vars=c("seqnames", "start", "end", "strand", 
                            "spliceSiteID"))
    
    # merge junction information with splice site annotation (theta)
    junction_dt <- merge(junction_dt, ss_dt[,.(spliceSiteID, sampleID, n)], 
                         all.x=TRUE, by.x=c("startID", "sampleID"), 
                         by.y=c("spliceSiteID", "sampleID"), sort=FALSE)
    setnames(junction_dt, "n", "n_donor")
    junction_dt <- merge(junction_dt, ss_dt[,.(spliceSiteID, sampleID, n)], 
                         all.x=TRUE, by.x=c("endID", "sampleID"), 
                         by.y=c("spliceSiteID", "sampleID"), sort=FALSE)
    setnames(junction_dt, "n", "n_acceptor")
    rm(ss_dt)
    gc()
    
    # deal with missing endIDs
    junction_dt_nacceptor <- data.table(
        as.data.table(rowRanges(fds, type="j"))[,.(startID, endID)], 
        as.matrix(N(fds, type="psi3")))
    junction_dt_nacceptor <- melt(junction_dt_nacceptor, 
                                  variable.name="sampleID", value.name="n_psi3", 
                                  id.vars=c("startID", "endID"))
    junction_dt[, n_psi3:=junction_dt_nacceptor[,n_psi3]]
    rm(junction_dt_nacceptor)
    gc()
    junction_dt_ndonor <- data.table(
        as.data.table(rowRanges(fds, type="j"))[,.(startID, endID)], 
        as.matrix(N(fds, type="psi5")))
    junction_dt_ndonor <- melt(junction_dt_ndonor, 
                               variable.name="sampleID", value.name="n_psi5", 
                               id.vars=c("startID", "endID"))
    junction_dt[, n_psi5:=junction_dt_ndonor[,n_psi5]]
    rm(junction_dt_ndonor)
    gc()
    
    # replace n (with non-split counts) by n_psi3/n_psi5 if NA
    junction_dt[is.na(n_acceptor), n_acceptor:=n_psi3]
    junction_dt[is.na(n_donor), n_donor:=n_psi5]
    junction_dt[n_acceptor < n_psi3, n_acceptor:=n_psi3]
    junction_dt[n_donor < n_psi5, n_donor:=n_psi5]
    
    # calculate phi
    junction_dt[, phi:= k / ((n_donor + n_acceptor)/2)]
    junction_dt[is.nan(phi), phi:=0] # n_donor = n_acceptor = 0
    
    # convert to matrix to store it as assay in the fds
    phi <- matrix(junction_dt[,phi], nrow=nrow(fds), ncol=ncol(fds), byrow=FALSE)
    rownames(phi) <- rownames(fds)
    colnames(phi) <- colnames(fds)
    assay(fds, "phi", type="psi5", withDimnames=FALSE) <- phi
    return(fds)
}
