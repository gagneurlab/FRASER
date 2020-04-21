context("testing getter/setter functions")

test_that("counts", {
    fds <- getFraser()
    expect_equal(counts(fds, type="psi5"),    K(fds, "psi5"))
    expect_equal(counts(fds, type="psi3"),    K(fds, "psi3"))
    expect_equal(counts(fds, type="psiSite"), K(fds, "psiSite"))
    expect_equal(as.matrix(counts(fds, type="psi5", side='other')),
            as.matrix(N(fds, "psi5") - K(fds, "psi5")))
    expect_equal(as.matrix(counts(fds, type="psi3", side='other')),
            as.matrix(N(fds, "psi3") - K(fds, "psi3")))
    expect_equal(as.matrix(counts(fds, type="psiSite", side='other')),
            as.matrix(N(fds, "psiSite") - K(fds, "psiSite")))
})
