context("Test stats calculations")

test_that("PSI value calculation", {
    
    fds <- getFraser()
    
    psiVal <- K(fds, "psi5") / N(fds, "psi5")
    
    # no pseudo counts in psi values
    expect_equal(psiVal[!is.na(psiVal)], assay(fds, "psi5")[!is.na(psiVal)])
    
    # NAs turned into 1
    expect_equal(assay(fds, "psi5")[is.na(psiVal)], rep(1, sum(is.na(psiVal))))
    
    # NAs where N==0
    expect_true(all(N(fds, "psi5")[!is.na(psiVal)] != 0))
    expect_true(all(N(fds, "psi5")[ is.na(psiVal)] == 0))
})

test_that("Zscore calculation", {
    fds <- getFraser(clean = TRUE)
    
    # prepare zScore input for logit scale
    psiVal <- (K(fds, "psi5") + pseudocount())/(N(fds, "psi5") + 2*pseudocount())
    mu <- predictedMeans(fds, "psi5")
    residual <- qlogis(psiVal) - qlogis(mu)
    
    # compute zscore
    zscores <- (residual - rowMeans(residual)) / rowSds(residual)
    
    expect_equal(zscores, zScores(fds, "psi5"))
})
