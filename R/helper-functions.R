########
## @author Christian Mertes \email{mertes@@in.tum.de}
##
## Helperfunctions to convert or extract data within the FraseR package
##

#'
#' clear the files in the cache to start fresh
#' @export
cleanCache <- function(fds, ...){
    stopifnot(class(fds) == "FraseRDataSet")
    # clean cache
    cacheDir <- file.path(workingDir(fds), "cache")
    if(dir.exists(cacheDir)){
        unlink(cacheDir, recursive=TRUE)
    }
}

#'
#' checks if the given type is part of the correct category
#' to map correctly to a read type
#'
#' @noRd
checkReadType <- function(fds, type){
    # check if type is null or missing
    if(missing(type) | is.null(type)){
        warning("Read type was not specified! We will assume the default: 'j'")
        return("j")
    }

    stopifnot(isScalarCharacter(type))
    correctTypes <- c(psi3="j", psi5="j", psiSite="ss")

    # check if it is already the correct type
    if(type %in% correctTypes) return(type)

    # check if psitype is given
    if(type %in% names(correctTypes)) return(correctTypes[type])

    # check assay names
    atype <- whichReadType(fds, type)
    if(isCorrectType(atype)) return(atype)

    stop("Given read type: '", type, "' not recognized. ",
            "It needs to be 'j' (junction) or 'ss' (splice sites)",
            "\nor an existing assay name within the given object."
    )
}

#'
#' returns the read type based on the given assay name
#'
#' @noRd
whichReadType <- function(fds, name){
    stopifnot(isScalarCharacter(name))
    fdsNames <- assayNames(fds)
    if(!name %in% fdsNames){
        return(NA)
    }
    nsrNamesL <- length(assayNames(nonSplicedReads(fds)))
    fdsNamesL <- length(fdsNames)

    return(ifelse(
        which(fdsNames == name) <= fdsNamesL - nsrNamesL,
        "j",
        "ss"
    ))
}

#'
#' Removes the white spaces to have a cleaner file path
#'@noRd
nameNoSpace <- function(name){
    if(class(name) == "FraseRDataSet") name <- name(name)
    stopifnot(isScalarCharacter(name))
    gsub("\\s+", "_", name, perl=TRUE)
}

#'
#' convert the input of NA to FALSE
#' @noRd
na2false <- function(x){
    if(any(class(x) %in% c("DataFrame", "matrix", "data.frame"))){
        stopifnot(dim(x)[2] == 1)
        x <- as.vector(as.matrix(x)[,1])
    }
    x[is.na(x)] <- FALSE
    return(x)
}

#'
#' the qq plot function
#' @noRd
fraserQQplot <- function(p, main=NULL, ...){
    zeroOffset <- 10e-100

    # my observerd and expected values
    o = -log10(sort(p + zeroOffset, decreasing=FALSE))
    e = -log10(ppoints(length(p)))[1:length(o)]

    # plot it
    plot(e, o, pch=19, cex=1, main=main, ...,
         xlab=expression(Expected~~-log[10](italic(p))),
         ylab=expression(Observed~~-log[10](italic(p))),
         xlim=c(0,max(e)), ylim=c(0,max(o))
    )
    lines(e, e, col="red")
}

#'
#' TODO this is not used yet or documented
#'
#' @noRd
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



