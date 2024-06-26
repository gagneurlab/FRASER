context("Test FRASER pipeline")

test_that("FRASER function", {
    fds <- createTestFraserDataSet()
    expect_is(fds, "FraserDataSet")
    anames <- c(psiTypes, paste0(c("delta", "predictedMeans", 
            "pvaluesBetaBinomial", "padjBetaBinomial"), "_", 
            rep(psiTypes, 5)))
    expect_equal(anames %in% assayNames(fds), !logical(length(anames)))
})


