#
# all functions to generate test objects
#
test_generate_count_example <- function(recount=FALSE){
    # get a new object for only one sample
    if(recount || !"test_fdsSample3" %in% ls()){

        test_fdsSample3 <- getFraser()[,"sample3",,drop=FALSE]
        name(test_fdsSample3) <- "onlySample3"

        # count the sample
        test_fdsSample3 <- countRNAData(test_fdsSample3, filter=FALSE, 
                                        recount=recount)
        test_fdsSample3 <- calculatePSIValues(test_fdsSample3)
    }

    #
    # create test objects
    #
    test_range <- GRanges(seqnames = "chr19", ranges = IRanges(
        start=c(7592515, 7594599, 7594599),
        end  =c(7592749, 7595171, 7595320)
    ))
    test_rangeOV  <- findOverlaps(test_range, test_fdsSample3, type = "equal")
    test_rangeFDS <- test_fdsSample3[to(test_rangeOV)]

    #
    # This is manually counted from the IGV browser
    #
    test_rawCountsJ   <- c(3, 13, 1)
    test_rawCountsSS  <- c(9, 10, 0, 0, 10)
    test_p5rawOCounts <- c(0, 1, 13)
    test_p3rawOCounts <- c(0, 0, 0)
    test_pSrawOCounts <- c(3, 3, 14, 13, 1)

    return(list(
        test_fdsSample3=test_fdsSample3,
        test_range=test_range,
        test_rangeOV=test_rangeOV,
        test_rangeFDS=test_rangeFDS,
        test_rawCountsJ=test_rawCountsJ,
        test_rawCountsSS=test_rawCountsSS,
        test_p5rawOCounts=test_p5rawOCounts,
        test_p3rawOCounts=test_p3rawOCounts,
        test_pSrawOCounts=test_pSrawOCounts
    ))
}

test_generate_strand_specific_count_example <- function(recount=FALSE){
    # get a new object for only one sample
    if(recount || !"test_fdsSample3_stranded" %in% ls()){
        
        test_fdsSample3_stranded <- createTestFraserSettings(
            workingDir=file.path(tempdir(), "strandSpecific"))[,"sample3",,drop=FALSE]
        name(test_fdsSample3_stranded) <- "onlySample3_stranded"
        strandSpecific(test_fdsSample3_stranded) <- TRUE
        pairedEnd(test_fdsSample3_stranded) <- TRUE
        
        # count the sample
        test_fdsSample3_stranded <- countRNAData(test_fdsSample3_stranded, 
                                                filter=FALSE, recount=recount)
        test_fdsSample3_stranded <- calculatePSIValues(test_fdsSample3_stranded)
    }
    
    #
    # create test objects
    #
    test_range_stranded <- GRanges(seqnames = "chr19", ranges = IRanges(
        start=c(7592515, 7592515, 7594599, 7594599, 7594599, 7594599),
        end  =c(7592749, 7592749, 7595171, 7595171, 7595320, 7595320)
    ),  strand=c("+", "-", "+", "-", "+", "-"))
    test_rangeOV_stranded  <- findOverlaps(test_range_stranded, 
                                    test_fdsSample3_stranded, type = "equal")
    test_rangeFDS_stranded <- 
        test_fdsSample3_stranded[to(test_rangeOV_stranded)]
    
    #
    # This is manually counted from the IGV browser 
    # colored by first of pair strand
    #
    
    # currently:
    # split counts: reads pairs with both ends mapped counted
    # nonSplit: read pairs counted, not requiring both ends mapped
    
    # split counts and ranges:
    # GRanges object with 5 ranges and 2 metadata columns:
    #     seqnames          ranges strand |   startID     endID
    # [1]    chr19 7592515-7592749      + |         4        51
    # [2]    chr19 7592515-7592749      - |        17        64
    # [3]    chr19 7594599-7595171      + |        11        57
    # [4]    chr19 7594599-7595171      - |        23        70
    # [5]    chr19 7594599-7595320      + |        11        58
    test_rawCountsJ_stranded   <- c(1, 1, 2, 8, 1)
    # if individual reads (with mapped mate in bam) are counted: 
    # c(2, 1, 2, 10, 1)
    
    # nonSplit ranges and counts:
    # GRanges object with 9 ranges and 2 metadata columns:
    #     seqnames          ranges strand | spliceSiteID     type
    # [1]    chr19 7592514-7592515      + |            4    Donor
    # [2]    chr19 7592749-7592750      + |           51 Acceptor
    # [3]    chr19 7594598-7594599      + |           11    Donor
    # [4]    chr19 7595171-7595172      + |           57 Acceptor
    # [5]    chr19 7595320-7595321      + |           58 Acceptor
    # [6]    chr19 7592514-7592515      - |           17 Acceptor
    # [7]    chr19 7592749-7592750      - |           64    Donor
    # [8]    chr19 7594598-7594599      - |           23 Acceptor
    # [9]    chr19 7595171-7595172      - |           70    Donor
    # manual counted based on read pairs:
    test_rawCountsSS_stranded  <- c(1, 3, 0, 0, 2, 6, 5, 0, 0) 
    
    test_p5rawOCounts_stranded <- c(0, 0, 1, 0, 2)
    test_p3rawOCounts_stranded <- c(0, 0, 0, 0, 0)
    test_pSrawOCounts_stranded <- c(2, 2, 3, 2, 1, 1, 1, 10, 10)
    
    
    
    return(list(
        test_fdsSample3_stranded=test_fdsSample3_stranded,
        test_range_stranded=test_range_stranded,
        test_rangeOV_stranded=test_rangeOV_stranded,
        test_rangeFDS_stranded=test_rangeFDS_stranded,
        test_rawCountsJ_stranded=test_rawCountsJ_stranded,
        test_rawCountsSS_stranded=test_rawCountsSS_stranded,
        test_p5rawOCounts_stranded=test_p5rawOCounts_stranded,
        test_p3rawOCounts_stranded=test_p3rawOCounts_stranded,
        test_pSrawOCounts_stranded=test_pSrawOCounts_stranded
    ))
}
