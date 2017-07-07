context("Test FraseR pipeline")

test_that("FraseR function", {
    fds <- getFraseR(TRUE)
    fds <- FraseR(settings=fds)
    
    expect_is(fds, "FraseRDataSet")
    expect_equal(dim(fds), c(94, 12))
    expect_equal(dim(nonSplicedReads(fds)), c(111, 12))
})

test_that("FraseRDataSet create settings", {
    settings <- createTestFraseRDataSet()
    expect_is(settings, "FraseRDataSet")
})


