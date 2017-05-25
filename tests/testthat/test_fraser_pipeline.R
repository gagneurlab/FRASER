context("Test FraseR pipeline")

test_that("FraseRSettings", {
    expect_equal("ok", "ok")
    settings <- createTestFraseRDataSet()
    expect_is(settings, "FraseRDataSet")
})

test_that("FraseR", {
    fds <- getFraseR()
    fds <- FraseR(settings=fds)

    expect_is(fds, "FraseRDataSet")
    expect_equal(dim(fds), c(94, 12))
    expect_equal(dim(nonSplicedReads(fds)), c(111, 12))
})
