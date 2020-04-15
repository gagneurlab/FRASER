context("Test FraseR pipeline")

test_that("FraseR function", {
    # fds <- getFraseR()
    # fds1 <- FraseR(fds=fds, q=2, iterations=2)
    # 
    # expect_is(fds, "FraseRDataSet")
    # expect_equal(dim(pVals(fds1, "psi5")),    dim(pVals(fds, "psi5")))
    # expect_equal(dim(pVals(fds1, "psiSite")), dim(pVals(fds, "psiSite")))
})

test_that("FraseRDataSet create settings", {
    fds <- createTestFraseRDataSet()
    expect_is(fds, "FraseRDataSet")
    anames <- c(psiTypes, paste0(c("delta", "predictedMeans", 
            "pvaluesBetaBinomial", "pajdBetaBinomial", "zScores"), "_", 
            rep(psiTypes, 5)))
    expect_equal(anames %in% assayNames(fds), !logical(length(anames)))
})


