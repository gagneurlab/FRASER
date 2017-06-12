context("Test generation of results")

test_that("Main plotting function", {
    # get subset to speed up test
    fds <- getFraseR()

    resFile <- plotSampleResults(fds, sampleID="sample1", browseIt=FALSE)

    expect_true(file.exists(resFile))
})

test_that("Multi Sample plotting", {
    # get subset to speed up test
    fds <- getFraseR()

    resFiles <- plotSampleResults(fds, browseIt=FALSE)

    expect_true(all(sapply(resFiles, file.exists)))

    expect_equal(length(resFiles), sum(!is.na(condition(fds))))
    expect_equal(length(resFiles), 4)
})


