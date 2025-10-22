context("Annotations methods")

test_that("annotateRanges", {

    fds <- createTestFraserSettings()
    fds <- countRNAData(fds)

    # compute stats
    fds <- calculatePSIValues(fds)
    # filter junctions with low expression
    fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=20,
					  minDeltaPsi=0.0, filter=TRUE)

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    orgDb <- org.Hs.eg.db

    # remove a splice site
    subset_fds <- fds[-c(1), , by = "ss"]

    # check to see if it is correctly removed
    expect_false(1 %in% mcols(subset_fds, type="theta")$spliceSiteID)

    # annoate genes
    subset_fds <- annotateRangesWithTxDb(subset_fds, txdb=txdb, orgDb=orgDb)

    # check if the junction is annotated
    j_mcols <- mcols(subset_fds, type = "j")
    expect_false(is.na(j_mcols$hgnc_symbol[j_mcols$startID == 1]))
})
