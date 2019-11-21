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
#' @param overwriteCts FALSE or TRUE (the default) the total counts (aka N) will
#'              be recalculated based on the existing junction counts (aka K)
#' @return FraseRDataSet
#' @export
#' @examples
#'   fds <- countRNAData(createTestFraseRSettings())
#'   fds <- calculatePSIValues(fds)
calculatePSIValues <- function(fds, types=psiTypes, overwriteCts=FALSE, 
                    BPPARAM=bpparam()){
    # check input
    stopifnot(is(fds, "FraseRDataSet"))

    # calculate PSI value for each sample
    for(psiType in unique(sapply(types, whichReadType, fds=fds))){
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

            # return only the important fields
            countData[,.(o5,o3,psi3,psi5)]
        }
    )
    names(psiValues) <- samples(fds)

    # merge it to a DataFrame and assign it to our object
    assay(fds, type="j", "psi5") <- sapply(psiValues, "[[", "psi5")
    assay(fds, type="j", "psi3") <- sapply(psiValues, "[[", "psi3")

    if(isTRUE(overwriteCts)){
        assay(fds, type="j", "rawOtherCounts_psi5") <-
            sapply(psiValues, "[[", "o5")
        assay(fds, type="j", "rawOtherCounts_psi3") <-
            sapply(psiValues, "[[", "o3")
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

            return(list(
                so=sdata[,os],
                psiSite=sdata[,psiValue]
            ))
        }
    )
    names(psiSiteValues) <- samples(fds)

    # merge it to a DataFrame and assign it to our object
    assay(fds, type="ss", psiName) <- sapply(psiSiteValues, "[[", "psiSite")
    if(isTRUE(overwriteCts)){
        assay(fds, type="ss", psiROCName) <- sapply(psiSiteValues, "[[", "so")
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
