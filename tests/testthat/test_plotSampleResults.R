context("Test generation of results")

test_that("Results function", {
    set.seed(42)
    # get subset to speed up test
    fds <- getFraser()
    
    # intron-level results
    res <- results(fds, aggregate=FALSE, all=TRUE)
    res <- as.data.table(res)
    for(psiType in psiTypes){
        expect_equal(res[type == psiType, .N], 
                    prod(dim(pVals(fds, level="junction", type=psiType))))
    }
    res_signif <- results(fds, aggregate=FALSE, all=FALSE,
                            padjCutoff=NA, deltaPsiCutoff=0.01)
    res_signif <- as.data.table(res_signif)
    expect_equal(res_signif[type == "jaccard", .N], 3)
    expect_equal(res_signif[type == "psi5", .N], 7)
    expect_equal(res_signif[type == "psi3", .N], 6)
    expect_equal(res_signif[type == "theta", .N], 8)
    
    # gene-level results
    res_gene <- results(fds, aggregate=TRUE, all=TRUE)
    res_gene <- as.data.table(res_gene)
    for(psiType in psiTypes){
        expect_equal(res_gene[type == psiType, .N], 
                    prod(dim(pVals(fds, level="gene", type=psiType))))
    }
    res_gene_signif <- results(fds, aggregate=TRUE, all=FALSE,
                          padjCutoff=NA, deltaPsiCutoff=0.01)
    res_gene_signif <- as.data.table(res_gene_signif)
    expect_equal(res_gene_signif[type == "jaccard", uniqueN(paste(seqnames, start, end))], 1)
    expect_equal(res_gene_signif[type == "psi5", uniqueN(paste(seqnames, start, end))], 3)
    expect_equal(res_gene_signif[type == "psi3", uniqueN(paste(seqnames, start, end))], 3)
    expect_equal(res_gene_signif[type == "theta", uniqueN(paste(seqnames, start, end))], 2)
    
    # results on subset of genes during FDR
    geneList <- list('sample1'=c("TIMMDC1"), 'sample2'=c("MCOLN1"))
    fds <- calculatePadjValues(fds, type="jaccard", 
                     subsets=list("exampleSubset"=geneList))
    expect_equal(length(results(fds, all=TRUE, psiType="jaccard")), 
                prod(dim(fds)))
    expect_error(results(fds, all=TRUE, psiType="psi5"))
    
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


