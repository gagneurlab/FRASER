context("Test betabinom pvalue calculations")

test_that("Check test randomization is correct", {
    # get subset to speed up test
    fds <- getFraseR()
    fds <- fds[which(mcols(fds, type="psi3")$psi3_tested)[1:20]]
    name(fds) <- "betabinomTest"
    fds <- saveFraseRDataSet(fds)
    parallel(fds) <- MulticoreParam(5)

    for(psiType in c("psi3", "psi5", "psiSite")){
        pvalls <- sapply(1:2, function(i){
            aname <- paste0("pvalue_", psiType)
            tmpfds <- pvalueByBetaBinomialPerType(fds, aname, psiType, betabinVglmTest)
            as.matrix(as.data.table(assays(tmpfds)[[aname]]))
        })

        expect_equal(pvalls[[1]], pvalls[[1]])
    }

})
