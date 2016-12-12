##
## @author Christian Mertes \email{mertes@@in.tum.de}
## 
## This file contains all functions for reading in data
## especially aligned RNA sequencing data 
##

#' 
#' This method extracts and counts the split reads from
#' a RNA bam file
#' @param settings A FraseRSetting object with all the information 
#'                 how and what to count
#' @param internBPPARAM A BiocParallel param object to configure the 
#'                      parallel backend of the internal loop
#' @return FraseRDataSet 
#' @export
#' @examples
#'   counRNAData(createTestFraseRSettings())
countRNAData <- function(settings, internBPPARAM=SerialParam()){
   
	# Check input TODO
    stopifnot(class(settings) == "FraseRSettings")
    stopifnot(!is(internBPParam, "BiocParallelParam"))

    # count splitreads first
    message("Start counting the split reads ...")
    countList <- bplapply(settings@sampleData[,bamFiles], 
                          FUN=.countSplitReads, 
                          settings=settings,
                          BPPARAM=settings@parallel,
                          internBPPARAM=internBPPARAM
    )
    names(countList) <- settings@sampleData[,sampleIDs]
    counts <- .mergeCounts(countList, settings@parallel)

    # count the retained reads
    message("Start counting the non spliced reads ...")
    countList <- bplapply(settings@sampleData[,bamFiles], 
                          FUN=.countNonSplicedReads, 
                          settings=settings,
                          targets=granges(counts),
                          BPPARAM=settings@parallel,
                          internBPPARAM=internBPPARAM
    )
    names(countList) <- settings@sampleData[,sampleIDs]
    site_counts <- .mergeCounts(countList, settings@parallel)
    mcols(site_counts)$type=factor(countList[[1]]$type, 
                    levels = c("Acceptor", "Donor")
    )
    
    # return it
    return(list(
        splitReads=counts, 
        nonSplicedReads=site_counts
    ))
}

#'
#' extracts the chromosomes within the given bamFile
#' 
.extractChromosomes <- function(bamFile){
    names(scanBamHeader(path(bamFile))[[path(bamFile)]]$target)
}

#'
#' count all split reads in a bam file
#'
.countSplitReads <- function(bamFile, strandSpecific, internBPPARAM){
    
    # parallelize over chromosomes
    chromosomes <- .extractChromosomes(bamFile)
    
    # extract the counts per chromosome
    countsList <- bplapply(chromosomes, FUN=.countSplitReadsPerChromosome,
                       bamFile=bamFile, 
                       strandSpecific=strandSpecific,
                       BPPARAM=internBPPARAM
    )
    
    # sort and merge the results befor returning
    return(sort(unlist(GRangesList(countsList))))
}

#'
#' counting the split reads per chromosome
#'
.countSplitReadsPerChromosome <- function(chromosome, bamFile, strandSpecific){
    # restrict to the chromosome only
    param <- ScanBamParam(which=GRanges(
            seqnames=chromosome, 
            ranges=IRanges(0, 536870912)
    ))

    # get reads from bam file
    galignment <- readGAlignments(bamFile, param=param)

    # remove the strand information if unstranded data
    if(!strandSpecific){
        strand(galignment) <- "*"
    }

    # dont count if there is nothing to count
    if(length(galignment) == 0){
        return(GRanges())	
    }

    # get the junction positions and their counts
    junctions <- unlist(junctions(galignment))
    junctionsCounts <- unique(junctions)

    # dont count anything if there is nothing to count
    if(length(junctionsCounts) == 0){
        return(junctionsCounts)
    }

    # count the data
    mcols(junctionsCounts)$count <- countOverlaps(junctionsCounts,
                                                  junctions,
                                                  type = 'equal'
    )

    # sort it and return the GRange object
    return(sort(junctionsCounts))
}
    

#' 
#' merge multiple samples into one SummarizedExperiment object
#' 
#' TODO: how to init non found junctions/splice sites (0L or NA)? 
#' 
.mergeCounts <- function(counts, BPPARAM=SerialParam()){
  
  # prepare range object
  sample_names <- names(counts)
  counts <- GRangesList(counts)
  ranges <- sort(unique(unlist(counts)))
  mcols(ranges)$count <- NULL
  names(ranges) <- NULL
  
  # merge each sample counts into the combined range object
  sample_counts <- bplapply(1:length(counts), ranges = ranges, 
                            counts = counts, BPPARAM=BPPARAM,
          FUN = function(i, ranges, counts){
                require(GenomicRanges)
            
                # get sample name
                sample_name <- names(counts)[i]
                
                ## TODO init with NA since we dont extracted
                # we are missing counts for junctions spliced in other samples
                sample_count <- rep(0L, length(ranges))
                
                # get overlap and add counts to the corresponding ranges
                overlaps <- findOverlaps(counts[[i]], ranges, type = "equal")
                sample_count[overlaps@to] <- mcols(counts[[i]])$count
                #mcols(ranges[overlaps@to,])[[sample_name]] <- mcols(counts[[i]])$count
                
                return(sample_count)
          }
  )
  
  # convert it to a DataFrame
  sample_counts_df <- DataFrame(
    matrix(unlist(sample_counts), ncol = length(sample_counts))
  )
  
  # merge it with the type columen and add it to the range object
  mcols(ranges) <- sample_counts_df
  
  # set correct naming 
  colnames(mcols(ranges)) <- sample_names
  
  # return the object
  return(ranges)
  
}

#' 
#' 
#' 
count_rna_data <- function(bamFile, strandSpecific=FALSE,
                           BPPARAM=SerialPAram()){
    
    # parallelize over chromosomes
	chromosomes <- names(scanBamHeader(bamfile)[[bamfile]]$target)
	
	# extract the counts per chromosome
	counts <- bplapply(chromosomes, bamFile=bamFile, 
                       strandSpecific=strandSpecific, BPPARAM=BPPARAM,
			FUN = function(chromosome, file, strand_specific){
				# restrict to the chromosome only
				param <- ScanBamParam(which = GRanges(
								seqnames = chromosome, 
								ranges = IRanges(0, 536870912)
						))
				
				# get reads from bam file
				galignment <- readGAlignments(bamfile, param = param)
				
				# if strand specific counting is not wanted remove the information
				if(!strand_specific){
					strand(galignment) <- "*"	
				}
				
				# dont count if there is nothing to count
				if(length(galignment) == 0){
					return(GRanges())	
				}
				
				# count non spliced reads @ junctions
				non_spliced_read_counts <- count_non_spliced_reads(galignment, junction_counts, strand_specific)
				
				# merge count data
				return(sort(unlist(
										GRangesList(junction_counts, non_spliced_read_counts)
								)))
			})
	
	# merge the chromosomes
	return(sort(unlist(GRangesList(counts))))
}



#'
#' counts non spliced reads based on the given target (acceptor/donor) regions
#'
.countNonSplicedReads <- function(bamFile, targets, strandSpecific=FALSE, internBPPARAM=SerialParam()){
    
    
    # extract donor and acceptor sites
    splice_site_coordinates <- extract_splice_site_coordinates(targets, strandSpecific)
    
    # extract the counts per chromosome
    countsList <- bplapply(splice_site_coordinates, bamFile=bamFile, 
                           strandSpecific=strandSpecific,
                           BPPARAM=internBPPARAM,
                           FUN=function(range, bamFile, strandSpecific){
                               single_read_fragments <- readGAlignments(bamFile, 
                                                param=ScanBamParam(which = range)) %>% 
                                                grglist() %>% reduce()
                               
                               countOverlaps(range, single_read_fragments, minoverlap = 2)
                           }
    )
	mcols(splice_site_coordinates)$count <- unlist(countsList)
	
	return(sort(splice_site_coordinates))
}


#'
#' 
#' 
extract_splice_site_coordinates <- function(junctions_gr, strandSpecific=FALSE){
	if(strandSpecific){
		splice_site_coords <- unlist(GRangesList(
				extract_splice_site_coordinates_per_strand(junctions_gr, "+"),
				extract_splice_site_coordinates_per_strand(junctions_gr, "-")
		))
	} else { 
		strand(junctions_gr) <- "*"
		splice_site_coords <- extract_splice_site_coordinates_per_strand(junctions_gr, "*")
	}
	
	return(sort(unique(splice_site_coords)))
}


#'
#' 
#' 
extract_splice_site_coordinates_per_strand <- function(junctions_gr, strand){
	
	# get only the correct strand features
	junctions_gr <- junctions_gr[strand(junctions_gr) == strand,]
	
	# left side (acceptor on + and donor on -) 
	left_side <- GRanges(
			seqnames = seqnames(junctions_gr),
			strand = strand(junctions_gr),
			ranges = IRanges(
					start = start(junctions_gr) - 1,
					end   = start(junctions_gr)
			),
			seqlengths = seqlengths(junctions_gr),
			seqinfo = seqinfo(junctions_gr)
	)
	
	# right side (acceptor on - and donor on +)
	right_side <- GRanges(
			seqnames = seqnames(junctions_gr),
			strand = strand(junctions_gr),
			ranges = IRanges(
					start = end(junctions_gr),
					end   = end(junctions_gr) + 1
			),
			seqlengths = seqlengths(junctions_gr),
			seqinfo = seqinfo(junctions_gr)
	)
	
	# annotate donor and acceptor sites
	if(strand == "+" | strand == "*"){
		mcols(left_side)$type = "Donor"
		mcols(right_side)$type = "Acceptor"
	} else {
		mcols(left_side)$type = "Acceptor"
		mcols(right_side)$type = "Donor"
	}
	
	return(sort(unlist(GRangesList(left_side, right_side))))
}


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
					)[,.(counts = sum(count,na.rm = T), is.na = all(is.na(count))), by = from]
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
calculate_psi_values <- function(count_data, BPPARAM = MulticoreParam(6, progressbar=T)){
	
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
			data[type == "Junction",psi_val:=get(sample)/sum(get(sample), na.rm = T),
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
get_assay_as_data_table <- function(se, assay, na_as_zero = T){
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
