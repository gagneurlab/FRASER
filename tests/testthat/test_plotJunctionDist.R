context("Test distribution plots for given results/junction")

test_that("Main junction distribution plot", {
    # get results
    fds <- getFraseR()
    res <- results(fds)

    # plot distributions
    expect_silent(plotJunctionDistribution(fds, res[res$type == "psi5"][1]))
    expect_silent(plotJunctionDistribution(fds, res[res$type == "psi3"][1]))
    expect_silent(plotJunctionDistribution(fds, res[res$type == "psiSite"][1]))
})

