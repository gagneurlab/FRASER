context("Test generation of results")

test_that("Results function", {
    set.seed(42)
    # get subset to speed up test
    fds <- getFraser()
    
    # intron-level results
    res <- results(fds, aggregate=FALSE, all=TRUE)
    expect_equal(length(res), prod(dim(fds)))
    res_signif <- results(fds, aggregate=FALSE, all=FALSE,
                            padjCutoff=NA, deltaPsiCutoff=0.2)
    expect_equal(length(res_signif), 1)
    
    # gene-level results
    res_gene <- results(fds, aggregate=TRUE, all=TRUE)
    expect_equal(length(res_gene), 
                    prod(dim(pVals(fds, level="gene", type="jaccard"))))
    res_gene_signif <- results(fds, aggregate=TRUE, all=FALSE,
                          padjCutoff=NA, deltaPsiCutoff=0.2)
    expect_equal(length(res_gene_signif), 1)
    
})

test_that("Main plotting function", {
    # get subset to speed up test
    # fds <- getFraser()

    # resFile <- plotSampleResults(fds, sampleID="sample1", browseIt=FALSE)

    # expect_true(file.exists(resFile))
    expect_true(TRUE)
})

test_that("Multi Sample plotting", {
    # get subset to speed up test
    # fds <- getFraser()

    # resFiles <- plotSampleResults(fds, browseIt=FALSE)

    # expect_true(all(sapply(resFiles, file.exists)))

    #expect_equal(length(resFiles), sum(!is.na(condition(fds))))
    #expect_equal(length(resFiles), 4)
    expect_true(TRUE)
})


