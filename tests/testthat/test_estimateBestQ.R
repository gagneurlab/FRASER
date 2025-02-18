context("Test estimation of optimal latent space dimension")

test_that("Test hyper param optimization", {
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

test_that("Test Optimal Hard Thresholding", {
  # Test error and warnings for special cases
  fds <- makeSimulatedFraserDataSet(m=2, j=2, dist="BB")
  fds <- calculatePSIValues(fds)
  
  expect_error(estimateBestQ(fds, type="psi5", useOHT=TRUE), 
               paste("OHT is only supported for the Intron Jaccard Index.",
                     "Set useOHT=FALSE to perform a hyperparameter optimization",
                     "for other metrics.", collapse = "\n"))
  expect_warning(estimateBestQ(fds, type="jaccard", useOHT=TRUE),
                 paste("Optimal latent space dimension is smaller than 2. Check your count matrices and",
                       "verify that all samples have the expected number of counts",
                       "For now, the latent space dimension is set to 2.", collapse = "\n"))
  
  # Test output for a regular example
  fds <- makeSimulatedFraserDataSet(m=15, j=20, dist="BB")
  fds <- calculatePSIValues(fds)
  fds <- estimateBestQ(fds, type="jaccard", useOHT=TRUE)
  
  expect_equal(c(1,3), dim(hyperParams(fds, type="jaccard")))
  expect_equal(metadata(fds)$hyperParams_jaccard[1,q],
               bestQ(fds, type="jaccard"))
  expect_is(suppressWarnings(plotEncDimSearch(fds, "jaccard", plotType="sv")), "ggplot")
})

test_that("Test optimalSVHTCoef", {
  expect_equal(optimalSVHTCoef(0.5),
               1.98, tolerance = 0.01)
  expect_equal(optimalSVHTCoef(0.1),
               1.58, tolerance = 0.01)
})

test_that("Test medianMarchenkoPastur", {
  # Expected outputs are derived from Table IV in Gavish and Donoho (2014)
  expect_equal(optimalSVHTCoef(0.5) / sqrt(medianMarchenkoPastur(100, 200)),
               2.1711, tolerance = 0.0001)
  expect_equal(optimalSVHTCoef(0.1) / sqrt(medianMarchenkoPastur(100, 1000)),
               1.6089, tolerance = 0.0001)
})
