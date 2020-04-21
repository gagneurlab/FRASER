context("Test stats calculations")

test_that("PSI value calculation", {
    # get subset to speed up test
    #fds <- getFraser()[,1]
    #fds <- fds[mcols(fds, type="psi3")$startID %in% c(1,2,6,7),]
    
    # expect it from the toy data after subsetting
    #expect_equal(length(fds), 10)
    #expect_equal(length(nonSplicedReads(fds)), 12)

    # set predifiend values
    #assays(fds)[["rawCountsJ"]]  <- rep(c(0,5), c(1,9))
    #assays(fds)[["rawCountsSS"]] <- rep(c(0,5), c(1,11))

    # calculate psi values
    #fds <- calculatePSIValues(fds)

    #assays(fds)[["psiSite"]]
    #counts(fds, type="psi3")

    #expect_equal("ok", "ok")
    #settings <- createTestFraserDataSet()
    #expect_is(settings, "FraserDataSet")
})

test_that("Zscore calculation", {
    #fds <- getFraser()
    #fds <- FRASER(fds)

    #expect_is(fds, "FraserDataSet")
    #expect_equal(dim(fds), c(94, 12))
    #expect_equal(dim(nonSplicedReads(fds)), c(111, 12))
})
