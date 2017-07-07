context("testing getter/setter functions")

test_that("counts", {
    fds <- getFraseR()
    expect_is(counts(fds, type="psi3"),    "DelayedMatrix")
    expect_is(counts(fds, type="psi5"),    "DelayedMatrix")
    expect_is(counts(fds, type="psiSite"), "DelayedMatrix")
    expect_is(counts(fds, type="psi3", side='other'),    "DelayedMatrix")
    expect_is(counts(fds, type="psi5", side='other'),    "DelayedMatrix")
    expect_is(counts(fds, type="psiSite", side='other'), "DelayedMatrix")
})
