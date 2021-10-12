context("Test FRASER pipeline")

test_that("FRASER function", {
    fds <- createTestFraserDataSet()
    expect_is(fds, "FraserDataSet")
    anames <- c(psiTypes, paste0(c("delta", "predictedMeans", 
            "pvaluesBetaBinomial_rho0.1", "padjBetaBinomial_rho0.1", 
            "zScores"), "_", 
            rep(psiTypes, 5)))
    expect_equal(anames %in% assayNames(fds), !logical(length(anames)))
})


