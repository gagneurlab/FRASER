context("Test distribution plots for given results/junction")

test_that("Main junction distribution plot", {
    # get results
    fds <- getFraseR()
    res <- results(fds, padjCutoff=1, zScoreCutoff=NA, deltaPsiCutoff=0.1)

    # plot distributions
    expect_silent(plotExpression(fds, result=res[1]))
    expect_silent(plotVolcano(fds, "sample1", "psi5"))
    expect_silent(plotExpectedVsObservedPsi(fds, result=res[2]))
})

