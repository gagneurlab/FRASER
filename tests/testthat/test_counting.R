context("Test counting")

fdsSample3 <- getFraseR()[,"sample3"]
name(fdsSample3) <- "onlySample3"
fdsSample3 <- countRNAData(fdsSample3)
fdsSample3 <- calculatePSIValues(fdsSample3)
testRange <- GRanges(seqnames = "chr19", ranges = IRanges(
    start=c(7592515, 7594599, 7594599),
    end  =c(7592749, 7595171, 7595320)
))
testRangeOV <- findOverlaps(testRange, fdsSample3, type = "equal")
testRangeFDS <- fdsSample3[to(testRangeOV)]

#
# This is manually counted from the IGV browser
#
rawCountsJ <- c(3, 13, 1)
rawCountsSS <- c(9, 10, 0, 0, 10)
p3rawOCounts <- c(0, 1, 13)
p5rawOCounts <- c(0, 0, 0)
pSrawOCounts <- c(3, 3, 14, 13, 1)


test_that("Count junctions", {
   expect_is(fdsSample3, "FraseRDataSet")

    # test how many ranges we found
    expect_equal(length(testRangeFDS), 3)
    expect_equal(length(nonSplicedReads(testRangeFDS)), 5)

    # test the manually counted positions
    expect_equal(as.vector(counts(testRangeFDS, type="j")), rawCountsJ)
    expect_equal(as.vector(counts(testRangeFDS, type="ss")), rawCountsSS)
})

test_that("Test psi values", {

    expect_equal(as.vector(counts(testRangeFDS, type="psi3")), rawCountsJ)
    expect_equal(p3rawOCounts,
        as.vector(counts(testRangeFDS, type="psi3", side="other"))
    )

    expect_equal(as.vector(counts(testRangeFDS, type="psi5")), rawCountsJ)
    expect_equal(p5rawOCounts,
        as.vector(counts(testRangeFDS, type="psi5", side="other"))
    )

    expect_equal(as.vector(counts(testRangeFDS, type="psiSite")), rawCountsSS)
    expect_equal(pSrawOCounts,
        as.vector(counts(testRangeFDS, type="psiSite", side="other"))
    )

    expect_equal(as.vector(assays(testRangeFDS)[["psi3"]]),
        rawCountsJ / (rawCountsJ + p3rawOCounts)
    )

    expect_equal(as.vector(assays(testRangeFDS)[["psi5"]]),
        rawCountsJ / (rawCountsJ + p5rawOCounts)
    )

    expect_equal(as.vector(assays(testRangeFDS)[["psiSite"]]),
                 rawCountsSS / (rawCountsSS + pSrawOCounts)
    )
})
