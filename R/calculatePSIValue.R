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
#'   fds <- countRNAData(createTestFraseRSettings())
#'   fds <- calculatePSIValues(fds)
calculatePSIValues <- function(fds, overwriteCts=FALSE){
    # check input
    stopifnot(class(fds) == "FraseRDataSet")

    # calculate PSI value for each sample
    for(psiType in c("psi5", "psi3", "psiSite")){
        if(!assayExists(fds, psiType)){
            fds <- calculatePSIValuePrimeSite(fds, psiType=psiType,
                    overwriteCts=overwriteCts)
            fds <- saveFraseRDataSet(fds)
            gc()
        }
    }

    # calculate the delta psi value
    for(psiType in c("psi3", "psi5", "psiSite")){
        assayName <- paste0("delta_", psiType)
        if(!assayExists(fds, assayName)){
            fds <- calculateDeltaPsiValue(fds, psiType, assayName)
            fds <- saveFraseRDataSet(fds)
            gc()
        }
    }

    # return it
    return(fds)
}


#'
#' calculates the PSI value for the given prime site of the junction
#'
#' @noRd
#'
#' we can get psi much simpler:
#' counts(fds, type='psi3', side='ofInterest')/
#' ( counts(fds, type='psi3', side='ofInterest')+ counts(fds, type='psi3', side='otherSide'))
#'
#'
calculatePSIValuePrimeSite <- function(fds, psiType, overwriteCts, logit=FALSE){
    stopifnot(class(fds) == "FraseRDataSet")
    stopifnot(isScalarCharacter(psiType))
    stopifnot(psiType %in% c("psi5", "psi3", "psiSite"))

    if(psiType=="psiSite"){
        return(calculateSitePSIValue(fds, overwriteCts))
    }

    message(date(), ": Calculate the PSI", psiType, " values ...")

    # convert psi type to the position of interest
    psiCol <- ifelse(psiType == "psi5", "start", "end")
    psiROCName <- paste0("rawOtherCounts_", psiType)
    if(! psiROCName %in% assayNames(fds)){
        overwriteCts <- TRUE
    }

    # generate a data.table from granges
    countData <- data.table(
        chr = as.factor(seqnames(fds)),
        start = start(fds),
        end = end(fds),
        strand = as.factor(strand(fds)),
        counts = NA
    )

    # calculate psi value
    psiValues <- bplapply(samples(fds), countData=countData, psiCol=psiCol,
        overwriteCts=overwriteCts, psiType=psiType, BPPARAM=parallel(fds), FUN =
                    function(sample, countData, psiCol, overwriteCts, psiType){
            suppressPackageStartupMessages(library(FraseR))

            # add sample specific counts to the data.table
            sample <- as.character(sample)
            sampleCounts <- as(counts(fds, type="j")[,sample], 'matrix')
            countData[,c(sample):=list(sampleCounts)]

            # calculate other split read counts
            if(isFALSE(overwriteCts)){
                countData[,rawOtherCounts:=as.vector(counts(fds, type=psiType,
                        side="oth")[,sample])]
            } else {
                countData[,rawOtherCounts:=sum(get(sample)) - get(sample),
                        by=eval(paste0("chr,", psiCol, ",strand"))
                ]
            }

            # calculate psi value
            if(isTRUE(logit)){
                countData[,psiValue:=logit(get(sample)/(get(sample) + rawOtherCounts + 1))]
            } else {
                countData[,psiValue:=get(sample)/(get(sample) + rawOtherCounts)]
            }

            # if psi is NA this means there were no reads at all so set it to 0
            countData[is.na(psiValue),psiValue:=0]

            return(list(
                    rawOtherCounts=countData[,rawOtherCounts],
                    psiValue=countData[,psiValue]
            ))
        }
    )
    names(psiValues) <- samples(fds)

    # merge it to a DataFrame and assign it to our object
    assays(fds, type="j")[[psiType]] <- DataFrame(
        lapply(psiValues, "[[", "psiValue")
    )
    if(isTRUE(overwriteCts)){
        assays(fds, type="j")[[psiROCName]] <- DataFrame(
            lapply(psiValues, "[[", "rawOtherCounts")
        )
    }

    return(fds)
}


#'
#' This function calculates the site PSI values for each splice site
#' based on the FraseRDataSet object
#'
#' @noRd
calculateSitePSIValue <- function(fds, overwriteCts){

    # check input
    stopifnot(class(fds) == "FraseRDataSet")

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
        overwriteCts=overwriteCts, BPPARAM=parallel(fds), FUN =
                    function(sample, countData, fds, overwriteCts){

            # add sample specific counts to the data.table
            sample <- as.character(sample)
            sampleCounts <- as.data.table(c(
                rep(as.vector(assays(fds)[["rawCountsJ"]][,sample]), 2),
                as.vector(assays(fds)[["rawCountsSS"]][,sample])
            ))
            colnames(sampleCounts) <- sample
            countData <- cbind(countData, sampleCounts)

            # calculate other split read counts
            if(isTRUE(overwriteCts)){
                countData[,rawOtherCounts:=sum(get(sample)) - get(sample),
                        by="spliceSiteID"]
            } else {
                countData[,rawOtherCounts:=as.vector(
                        counts(fds, type="psiSite", side="oth")[,sample])]
            }

            # remove the junction part since we only want to calculate the
            # psi values for the splice site itself
            countData <- countData[type=="spliceSite"]

            # calculate psi value
            countData[,psiValue:=get(sample)/(get(sample) + rawOtherCounts)]

            # if psi is NA this means there were no reads at all so set it to 0
            countData[is.na(psiValue),psiValue:=0]

            return(list(
                rawOtherCounts=countData[,rawOtherCounts],
                psiValue=countData[,psiValue]
            ))
        }
    )
    names(psiSiteValues) <- samples(fds)

    # merge it to a DataFrame and assign it to our object
    assays(fds, type="ss")[[psiName]] <- DataFrame(
        lapply(psiSiteValues, "[[", "psiValue"))
    if(isTRUE(overwriteCts)){
        assays(fds, type="ss")[[psiROCName]] <- DataFrame(
            lapply(psiSiteValues, "[[", "rawOtherCounts"))
    }

    return(fds)
}

#'
#' calculates the delta psi value and stores it as an assay
#'
calculateDeltaPsiValue <- function(fds, psiType, assayName){

    message(date(), ": Calculate the delta for ", psiType, " values ...")

    # get psi values
    psiVal <- assays(fds)[[psiType]]

    # psi - median(psi)
    deltaPsi  <- psiVal - rowMedians(psiVal, na.rm = TRUE)

    # use as.matrix to rewrite it as a new hdf5 array
    assays(fds, type=psiType)[[assayName]] <- deltaPsi

    return(fds)
}
