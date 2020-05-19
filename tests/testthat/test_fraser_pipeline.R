context("Test FRASER pipeline")

test_that("FRASER function", {
    # fds <- getFraser()
    # fds1 <- FRASER(fds=fds, q=2, iterations=2)
    # 
    # expect_is(fds, "FraserDataSet")
    # expect_equal(dim(pVals(fds1, "psi5")),    dim(pVals(fds, "psi5")))
    # expect_equal(dim(pVals(fds1, "psiSite")), dim(pVals(fds, "psiSite")))
})

test_that("FraserDataSet create settings", {
    fds <- createTestFraserDataSet()
    expect_is(fds, "FraserDataSet")
    anames <- c(psiTypes, paste0(c("delta", "predictedMeans", 
            "pvaluesBetaBinomial", "padjBetaBinomial", "zScores"), "_", 
            rep(psiTypes, 5)))
    expect_equal(anames %in% assayNames(fds), !logical(length(anames)))
})


