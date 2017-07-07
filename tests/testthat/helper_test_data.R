#
# all functions to generate test objects
#
test_generate_count_example <- function(){
    # get a new object for only one sample
    if("test_fdsSample3" %nin% ls()){

        test_fdsSample3 <<- getFraseR()[,"sample3"]
        name(test_fdsSample3) <<- "onlySample3"

        # count the sample
        test_fdsSample3 <- countRNAData(test_fdsSample3)
        test_fdsSample3 <- calculatePSIValues(test_fdsSample3)
    }

    #
    # create test objects
    #
    test_range <- GRanges(seqnames = "chr19", ranges = IRanges(
        start=c(7592515, 7594599, 7594599),
        end  =c(7592749, 7595171, 7595320)
    ))
    test_rangeOV <- findOverlaps(test_range, test_fdsSample3, type = "equal")
    test_rangeFDS <- test_fdsSample3[to(test_rangeOV)]

    #
    # This is manually counted from the IGV browser
    #
    test_rawCountsJ <- c(3, 13, 1)
    test_rawCountsSS <- c(9, 10, 0, 0, 10)
    test_p3rawOCounts <- c(0, 1, 13)
    test_p5rawOCounts <- c(0, 0, 0)
    test_pSrawOCounts <- c(3, 3, 14, 13, 1)

    return(list(
        test_fdsSample3=test_fdsSample3,
        test_range=test_range,
        test_rangeOV=test_rangeOV,
        test_rangeFDS=test_rangeFDS,
        test_rawCountsJ=test_rawCountsJ,
        test_rawCountsSS=test_rawCountsSS,
        test_p3rawOCounts=test_p3rawOCounts,
        test_p5rawOCounts=test_p5rawOCounts,
        test_pSrawOCounts=test_pSrawOCounts
    ))
}
