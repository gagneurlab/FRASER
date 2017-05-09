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
calculatePSIValues <- function(fds){
    # check input
    stopifnot(class(fds) == "FraseRDataSet")

    # calculate 3/5' PSI for each sample
    fds <- calculatePSIValuePrimeSite(fds, psiType="3")
    fds <- calculatePSIValuePrimeSite(fds, psiType="5")

    # calculate siteSplice values
    fds <- calculateSitePSIValue(fds)

    # return it
    return(fds)
}


#'
#' calculates the PSI value for the given prime site of the junction
#'
#' @noRd
calculatePSIValuePrimeSite <- function(fds, psiType){
    stopifnot(class(fds) == "FraseRDataSet")
    stopifnot(isScalarCharacter(psiType))
    stopifnot(psiType %in% c("3", "5"))

    message(date(), ": Calculate the PSI", psiType, " values ...")

    # convert psi type to the position of interest
    psiCol <- ifelse(psiType == "3", "start", "end")
    psiName <- paste0("psi", psiType)
    psiROCName <- paste0("rawOtherCounts_psi", psiType)

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
        BPPARAM=parallel(fds), FUN = function(sample, countData, psiCol){
            suppressPackageStartupMessages(library(FraseR))

            # add sample specific counts to the data.table
            sample <- as.character(sample)
            sampleCounts <- as.data.table(assays(fds)[["rawCountsJ"]][,sample])
            colnames(sampleCounts) <- sample
            countData <- cbind(countData, sampleCounts)

            # calculate other split read counts
            countData[,rawOtherCounts:=sum(get(sample)) - get(sample),
                    by=eval(paste0("chr,", psiCol, ",strand"))
            ]

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
    names(psiValues) <- samples(fds)

    # merge it to a DataFrame and assign it to our object
    assays(fds, type="j")[[psiName]] <- DataFrame(
        lapply(psiValues, "[[", "psiValue")
    )
    assays(fds, type="j")[[psiROCName]] <- DataFrame(
        lapply(psiValues, "[[", "rawOtherCounts")
    )

    return(fds)
}


#'
#' This function calculates the site PSI values for each splice site
#' based on the FraseRDataSet object
#'
#' @noRd
calculateSitePSIValue <- function(fds){

    # check input
    stopifnot(class(fds) == "FraseRDataSet")

    message(date(), ": Calculate the PSI site values ...")

    psiName <- "psiSite"
    psiROCName <- "rawOtherCounts_psiSite"

    # prepare data table for calculating the psi value
    countData <- data.table(
        spliceSiteID=c(
            rowData(fds)[["startID"]],
            rowData(fds)[["endID"]],
            rowData(nonSplicedReads(fds))[["spliceSiteID"]]
        ),
        type=rep(
            c("junction", "spliceSite"),
            c(length(fds)*2, length(nonSplicedReads(fds)))
        )
    )

    psiSiteValues <- bplapply(samples(fds), countData=countData, fds=fds,
        BPPARAM=parallel(fds), FUN = function(sample, countData, fds){

            # add sample specific counts to the data.table
            sample <- as.character(sample)
            sampleCounts <- as.data.table(c(
                rep(as.vector(assays(fds)[["rawCountsJ"]][,sample]), 2),
                as.vector(assays(fds)[["rawCountsSS"]][,sample])
            ))
            colnames(sampleCounts) <- sample
            countData <- cbind(countData, sampleCounts)

            # calculate other split read counts
            countData[,rawOtherCounts:=sum(get(sample)) - get(sample),
                    by="spliceSiteID"
            ]

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
        lapply(psiSiteValues, "[[", "psiValue")
    )
    assays(fds, type="ss")[[psiROCName]] <- DataFrame(
        lapply(psiSiteValues, "[[", "rawOtherCounts")
    )

    return(fds)
}


