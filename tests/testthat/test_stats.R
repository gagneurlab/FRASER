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

# test_that("Zscore calculation", {
#     fds <- getFraser(clean = TRUE)
#     
#     # prepare zScore input for logit scale
#     psiVal <- (K(fds, "jaccard") + pseudocount())/(N(fds, "jaccard") + 2*pseudocount())
#     mu <- predictedMeans(fds, "jaccard")
#     residual <- qlogis(psiVal) - qlogis(mu)
#     
#     # compute zscore
#     zscores <- (residual - rowMeans(residual)) / rowSds(residual)
#     
#     expect_equal(zscores, zScores(fds, "jaccard"))
# })

test_that("Gene p value calculation with NAs", {
    fds <- getFraser()
    fds <- fds[15:24,]
    mcols(fds, type="j")$hgnc_symbol <- rep(c("geneA", "geneB", "geneC"), 
                                            times=c(3, 4, 3))
    mcols(fds, type="ss")$hgnc_symbol <- rep(c("geneA", "geneB", "geneC"), 
                                            times=c(4, 6, 4))
    
    # simulate junction with bad rho fit to create partly NAs
    rho <- rho(fds, type="jaccard")
    rho[c(1, 4:7)] <- 0.5
    rho(fds, type="jaccard") <- rho
  
    # calc p values
    fds <- calculatePadjValues(fds, type="jaccard", rhoCutoff=0.1)
    
    # check dimension of junction-, site- and gene-level pval matrices
    expect_equal(nrow(pVals(fds, type="jaccard", level="junction")), nrow(fds))
    expect_equal(nrow(pVals(fds, type="jaccard", level="site", 
                            filters=list(rho=0.1))), nrow(fds))
    expect_equal(nrow(pVals(fds, type="jaccard", level="gene", 
                            filters=list(rho=0.1))), 3)
    
    # check jaccard pvals are partly NAs
    expect_true(all(is.na(pVals(fds, type="jaccard", level="site", 
                            filters=list(rho=0.1))[4:7,])))
    expect_true(all(is.na(pVals(fds, type="jaccard", level="gene", 
                            filters=list(rho=0.1))["geneB",])))
    expect_true(all(is.na(padjVals(fds, type="jaccard", level="site", 
                            filters=list(rho=0.1))[4:7,])))
    expect_true(all(is.na(padjVals(fds, type="jaccard", level="gene", 
                            filters=list(rho=0.1))["geneB",])))
    
    # simulate junction with bad rho fit to create partly NAs
    rho <- rho(fds, type="jaccard")
    rho <- rep(0.5, length(rho))
    rho(fds, type="jaccard") <- rho
    fds <- calculatePadjValues(fds, type="jaccard", rhoCutoff=0.1)
    
    # check jaccard pvals are all NAs
    expect_true(all(is.na(pVals(fds, type="jaccard", level="site", 
                                filters=list(rho=0.1)))))
    expect_true(all(is.na(pVals(fds, type="jaccard", level="gene", 
                                filters=list(rho=0.1)))))
    expect_true(all(is.na(padjVals(fds, type="jaccard", level="site", 
                                filters=list(rho=0.1)))))
    expect_true(all(is.na(padjVals(fds, type="jaccard", level="gene", 
                                filters=list(rho=0.1)))))
})

test_that("FDR on subset of genes", {
    fds <- getFraser()
    mcols(fds, type="j")$hgnc_symbol <- 
        rep(c("geneA", "geneB", "geneC", "geneD", "geneE"), 
            times=c(3, 7, 5, 4, 8))
    
    # define gene subset per sample
    genes_per_sample <- list(
        "sample1" = c("geneE", "geneC", "geneA"),
        "sample2" = c("geneB"),
        "sample3" = c("geneA", "geneB", "geneC", "geneD")
    )
    
    subsetName <- "subset_test"
    fds <- calculatePadjValuesOnSubset(fds, genesToTest=genes_per_sample, 
                                       subsetName=subsetName, type="jaccard")
    subset_padj <- padjVals(fds, type="jaccard", subsetName=subsetName)
    expect_true(is(subset_padj, "matrix"))
    expect_true(nrow(subset_padj) == 27)
    expect_true(ncol(subset_padj) == 3)
    subset_padj_gene <- padjVals(fds, type="jaccard", level="gene", 
                                    subsetName=subsetName)
    expect_true(is(subset_padj_gene, "matrix"))
    expect_true(nrow(subset_padj_gene) == 5)
    expect_true(ncol(subset_padj_gene) == 3)
      
})
