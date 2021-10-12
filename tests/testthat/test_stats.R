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

test_that("Gene p value calculation with NAs", {
    fds <- getFraser()
    fds <- fds[15:24,]
    mcols(fds, type="j")$hgnc_symbol <- rep(c("geneA", "geneB", "geneC"), 
                                            times=c(3, 4, 3))
    mcols(fds, type="ss")$hgnc_symbol <- rep(c("geneA", "geneB", "geneC"), 
                                            times=c(4, 6, 4))
    
    # simulate junction with bad rho fit
    rho_5 <- rho(fds, type="psi5")
    rho_5[c(1, 4:7)] <- 0.5
    rho(fds, type="psi5") <- rho_5
    
    rho_3 <- rho(fds, type="psi3")
    rho_3 <- rep(0.5, length(rho_3))
    rho(fds, type="psi3") <- rho_3
  
    # calc p values
    fds <- calculatePadjValues(fds, type="psi5", rhoCutoff=0.1)
    fds <- calculatePadjValues(fds, type="psi3", rhoCutoff=0.1)
    
    # check psi5 pvals are partly NAs
    expect_equal(pVals(fds, type="psi5", level="site", 
                        filters=list(rho=0.1))[4:7,1],
                 as.double(rep(NA, 4)))
    expect_equal(pVals(fds, type="psi5", level="gene", 
                        filters=list(rho=0.1))[4:7,2],
                 as.double(rep(NA, 4)))
    expect_equal(padjVals(fds, type="psi5", level="site", 
                        filters=list(rho=0.1))[4:7,1],
                 as.double(rep(NA, 4)))
    expect_equal(padjVals(fds, type="psi5", level="gene", 
                        filters=list(rho=0.1))[4:7,2],
                 as.double(rep(NA, 4)))
    
    # check psi3 pvals are all NAs
    expect_equal(pVals(fds, type="psi3", level="site", 
                        filters=list(rho=0.1))[,1],
                 as.double(rep(NA, nrow(fds))))
    expect_equal(pVals(fds, type="psi3", level="gene", 
                        filters=list(rho=0.1))[,2],
                 as.double(rep(NA, nrow(fds))))
    expect_equal(padjVals(fds, type="psi3", level="site", 
                        filters=list(rho=0.1))[,1],
                 as.double(rep(NA, nrow(fds))))
    expect_equal(padjVals(fds, type="psi3", level="gene", 
                        filters=list(rho=0.1))[,2],
                 as.double(rep(NA, nrow(fds))))
})
