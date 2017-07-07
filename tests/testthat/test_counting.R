context("Test counting")

test_that("Count junctions", {
    attach(test_generate_count_example())

    expect_is(test_fdsSample3, "FraseRDataSet")

    # test how many ranges we found
    expect_equal(length(test_rangeFDS), 3)
    expect_equal(length(nonSplicedReads(test_rangeFDS)), 5)

    # test the manually counted positions
    expect_equal(as.vector(counts(test_rangeFDS, type="j")), test_rawCountsJ)
    expect_equal(as.vector(counts(test_rangeFDS, type="ss")), test_rawCountsSS)
})

test_that("Test psi values", {
    attach(test_generate_count_example())

    expect_equal(as.vector(counts(test_rangeFDS, type="psi3")), test_rawCountsJ)
    expect_equal(test_p3rawOCounts,
        as.vector(counts(test_rangeFDS, type="psi3", side="other"))
    )

    expect_equal(as.vector(counts(test_rangeFDS, type="psi5")), test_rawCountsJ)
    expect_equal(test_p5rawOCounts,
        as.vector(counts(test_rangeFDS, type="psi5", side="other"))
    )

    expect_equal(as.vector(counts(test_rangeFDS, type="psiSite")), test_rawCountsSS)
    expect_equal(test_pSrawOCounts,
        as.vector(counts(test_rangeFDS, type="psiSite", side="other"))
    )

    expect_equal(as.vector(assays(test_rangeFDS)[["psi3"]]),
        test_rawCountsJ / (test_rawCountsJ + test_p3rawOCounts)
    )

    expect_equal(as.vector(assays(test_rangeFDS)[["psi5"]]),
        test_rawCountsJ / (test_rawCountsJ + test_p5rawOCounts)
    )

    expect_equal(as.vector(assays(test_rangeFDS)[["psiSite"]]),
        test_rawCountsSS / (test_rawCountsSS + test_pSrawOCounts)
    )
})
