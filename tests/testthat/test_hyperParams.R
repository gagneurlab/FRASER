context("Test hyper param optimization")

test_that("Test hyper param testing", {
    fds <- makeSimulatedFraserDataSet(m=15, j=20, dist="BB")
    fds <- calculatePSIValues(fds)
    
    # test BB no hyper params and accessors
    fds <- estimateBestQ(fds, type="psi3", useOHT=FALSE, implementation="BB")
    expect_true(is.na(hyperParams(fds, type="psi3")[,q]))
    
    fds <- estimateBestQ(fds, type="psi5", useOHT=FALSE, minDeltaPsi=0.01, plot=FALSE,
            q_param = c(2,3), noise_param=c(1), iterations=2, BPPARAM=bpparam())
    
    expect_equal(c(2,5), dim(hyperParams(fds, type="psi5", all=TRUE)))
    expect_equal(c(1,5), dim(hyperParams(fds, type="psi5", all=FALSE)))
    expect_equal(metadata(fds)$hyperParams_psi5[order(-aroc)][1,q],
            bestQ(fds, type="psi5"))
    expect_equal(1, bestNoise(fds, type="psi5"))
    expect_is(suppressWarnings(plotEncDimSearch(fds, "psi5")), "ggplot")
})
