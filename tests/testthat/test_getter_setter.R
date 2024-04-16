context("testing getter/setter functions")

test_that("counts", {
    fds <- createTestFraserDataSet()
    expect_equal(counts(fds, type="psi5"),    K(fds, "psi5"))
    expect_equal(counts(fds, type="psi3"),    K(fds, "psi3"))
    expect_equal(counts(fds, type="theta"), K(fds, "theta"))
    expect_equal(counts(fds, type="jaccard"),    K(fds, "jaccard"))
    expect_equal(counts(fds, type="psi5", side='other'), 
            N(fds, "psi5") - K(fds, "psi5"))
    expect_equal(counts(fds, type="psi3", side='other'),
            N(fds, "psi3") - K(fds, "psi3"))
    expect_equal(counts(fds, type="theta", side='other'),
            N(fds, "theta") - K(fds, "theta"))
    expect_equal(counts(fds, type="jaccard", side='other'),
                 N(fds, "jaccard") - K(fds, "jaccard"))
})

test_that("generic functions", {
    # test that the correct functions are selected and called
    methDef <- selectMethod("[", c("FraserDataSet", "ANY", "ANY", drop="ANY"))
    expect_equal(slot(methDef, "defined")[["x"]], "FraserDataSet")
    methDef <- selectMethod("[", 
            c("FraserDataSet", "ANY", "ANY", drop="missing"))
    expect_equal(slot(methDef, "defined")[["x"]], "FraserDataSet")
})
