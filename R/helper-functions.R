########
## @author Christian Mertes \email{mertes@@in.tum.de}
## 
## Helperfunctions to convert or extract data within the FraseR package
##


#'
#' convert a data.table to a DataFrame and keep the colname names
#' 
#' @noRd
.asDataFrame <- function(dataframe, colname = colnames(dataframe)){
    dataframe <- DataFrame(dataframe)
    colnames(dataframe) <- colname
    return(dataframe)
}


#'
#' get the assay as data.table from a SummarizedExperiment object
#' 
#' @noRd
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



