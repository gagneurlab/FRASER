context("Test betabinom pvalue calculations")

test_that("Check test randomization is correct", {
    # get subset to speed up test
    fds <- getFraseR()
    fds <- fds[which(mcols(fds1, type="psi3")$psi3_tested)[1:20]]
    name(fds) <- "betabinomTest"
    parallel(fds) <- MulticoreParam(10)

    for(psiType in c("psi3", "psi5", "psiSite")){
        fds1 <- testPsiWithBetaBinomialPerType(fds, psiType, betabinVglmTest)
        pv1 <- as(assays(fds1)[[paste0("pvalue_", psiType)]], "matrix")
        fds2 <- testPsiWithBetaBinomialPerType(fds, psiType, betabinVglmTest)
        pv2 <- as(assays(fds2)[[paste0("pvalue_", psiType)]], "matrix")

        expect_equal(pv1[,1], pv2[,1])
    }

})
