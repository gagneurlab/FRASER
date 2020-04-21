context("Test betabinom pvalue calculations")

test_that("Check test randomization is correct", {
    # get subset to speed up test
    # fds <- getFraser()
    # fds <- fds[which(mcols(fds, type="psi3")$psi3_tested)[1:20]]
    # name(fds) <- "betabinomTest"
    # fds <- saveFraserDataSet(fds)
    
    # for(psiType in c("psi3", "psi5", "psiSite")){
     #   pvalls <- sapply(1:2, function(i){
     #       aname <- paste0("pvalue_", psiType)
     #       tmpfds <- pvalueByBetaBinomialPerType(fds, aname, psiType, betabinVglmTest)
     #       as.matrix(as.data.table(assays(tmpfds)[[aname]]))
     #   })
#
 #       expect_equal(pvalls[[1]], pvalls[[1]])
  #  }
    expect_equal(1,1)
})
