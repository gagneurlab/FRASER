

#'
#' 
#' 
calculate_intron_retaintion <- function(data){
    junctions <- data[rowData(data)$type == "Junction"]
    
    for(splice_type in c("Acceptor", "Donor")){
        modified_rows <- rowData(data)$type == splice_type
        splice_site <- data[modified_rows]
        splice_ranges <- rowRanges(splice_site)
        
        # shift for start/end overlap
        splice_ranges <- shift(splice_ranges, ifelse(splice_type == "Acceptor", 1, -1))
        
        # find overlap
        overlap <- findOverlaps(splice_ranges, junctions, type = 
                                    ifelse(splice_type == "Acceptor", "start", "end")
        )
        
        # sum up the junctions per site per sample
        junction_counts_per_site <- mclapply(colData(data)$sample, overlap = overlap, 
                                             junctions = junctions, mc.cores = 4, 
                                             FUN = function(sample, overlap, junctions){
                                                 dt <- data.table(
                                                     from = overlap@from,
                                                     count = as.numeric(assays(junctions[,sample])$counts[,1])[overlap@to]
                                                 )[,.(counts = sum(count,na.rm = TRUE), is.na = all(is.na(count))), by = from]
                                                 return(dt[,counts])
                                             }
        )
        
        # add junction counts
        splice_site <- merge_assay("junction_site_counts", splice_site, 
                                   junction_counts_per_site, unique(overlap@from))
        
        # add junction persent spliced in value (junction_psi) 
        splice_site <- merge_assay("junction_site_psi", splice_site,
                                   get_assay_as_data_table(splice_site, "junction_site_counts") / (
                                       get_assay_as_data_table(splice_site, "junction_site_counts") +
                                           get_assay_as_data_table(splice_site, "counts")
                                   )
        )
        
        
        # wrote it to the data
        assays(data)$junction_site_counts <- DataFrame
        assays(data[modified_rows])$junction_site_counts <- assays(splice_site)$junction_site_counts
        assays(data[modified_rows])$junction_site_psi    <- assays(splice_site)$junction_site_psi
        
        assays(data)<- NULL
    }
    
    
    
}


#'
#' 
#' 
calculate_psi_values <- function(count_data, BPPARAM = MulticoreParam(6, progressbar=TRUE)){
    
    # generate a data.table from granges
    tmp_data <- cbind(
        data.table(
            chr = as.factor(seqnames(count_data)),
            start = start(count_data),
            end = end(count_data),
            strand = as.factor(strand(count_data)),
            type = rowRanges(count_data)$type,
            counts = NA
        ),
        as.data.table(assays(count_data)$counts)
    )
    
    # calculate 3/5' psi for each sample
    assays(count_data)$psi3 <- calculate_x_prime_psi_values(data = tmp_data,
                                                            samples = colData(count_data)$sample, psi_type = "3'", BPPARAM = BPPARAM
    )
    assays(count_data)$psi5 <- calculate_x_prime_psi_values(data = tmp_data,
                                                            samples = colData(count_data)$sample, psi_type = "5'", BPPARAM = BPPARAM
    )
    
    return(count_data)
    
}


#'
#' 
#' 
calculate_x_prime_psi_values <- function(data, samples, psi_type = "3'",
                                         BPPARAM = MulticoreParam(6)){
    
    # convert psi type to the position of interest 
    if(psi_type == "3'"){
        psi_col = "start"
    } else {
        psi_col = "end"
    }
    
    # calculate psi value
    psi_values <- mclapply(samples, data = data, psi_col = psi_col,
                           mc.cores = bpworkers(BPPARAM), mc.preschedule = FALSE, 
                           FUN = function(sample, data, psi_col){
                               
                               # check name, du to conversion
                               if(grepl("^\\d+$", sample)){
                                   sample <- paste0("X", sample)
                               }
                               
                               # init psi
                               data[,psi_val:=as.numeric(NA)]
                               
                               # calculate psi
                               data[type == "Junction",psi_val:=get(sample)/sum(get(sample), na.rm = TRUE),
                                    by = eval(paste0("chr,", psi_col,",strand"))
                                    ]
                               
                               return(data[,psi_val])
                           }
    )
    
    # merge it and set the column names
    df = DataFrame(matrix(unlist(psi_values), ncol = length(samples)))
    names(df) <- samples 
    
    return(df)
    
}


#' 
#' 
#' 
convert_dataframe_columns_to_Rle <- function(data, index2convert = 1:dim(dataframe)[2]){
    
    # convert all given indices
    for (i in index2convert){
        data[,i] <- Rle(data[,i])
    }
    
    return(data)
}


#'
#' get the assay as data.table
#' 
get_assay_as_data_table <- function(se, assay, na_as_zero = TRUE){
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
#' convert a data.table to a DataFrame and keep the colnames
#' 
as_DataFrame <- function(dataframe, colname = colnames(dataframe)){
    dataframe <- DataFrame(dataframe)
    colnames(dataframe) <- colnames
    return(dataframe)
}


#'
#' 
#' 
merge_assay <- function(assay1, se1, df1, rowidx1 = 1:dim(se1)[1]){
    # create empty dataframe for assay if not present yet
    if(!any(names(assays(se1)) %in% assay1)){
        tmp_df <- DataFrame(matrix(as.numeric(NA), nrow = dim(se1)[1], ncol = dim(se1)[2]))
        colnames(tmp_df) <- colData(se1)$sample
        assays(se1)[[assay1]] <- tmp_df
    } 
    
    # make a data frame out of the given data if needed
    if(!any(class(df1) %in% "DataFrame")){
        df1 <- as_DataFrame(df1, colData(se1)$sample)
    }
    
    # merge the data
    tmp_df <- assays(se1)[[assay1]]
    tmp_df[rowidx1,] <- df1
    assays(se1)[[assay1]] <- tmp_df
    
    # return it
    return(se1)
}

#'
#'
#' calculate the zscore for each psi value
calculate_zscore_splicing <- function(data, psi_type){
    
    # get raw data and replace NA's with zeros
    psi_val <- get_assay_as_data_table(data, psi_type)
    
    # z = ( x - mean ) / sd
    rowmean <- rowMeans(psi_val, na.rm = TRUE)
    rowsd   <- apply(psi_val, 1, sd, na.rm = TRUE)
    zscores = (psi_val - rowmean) / rowsd
    
    # add it to the SE object
    assays(data)[[paste0("zscore_", psi_type)]] <- 
        as_DataFrame(zscores, colData(data)$sample)
    
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
