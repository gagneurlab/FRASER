context("testing util functions")

test_that("get URL links", {
    dt <- data.table(
        seqnames = "chr1",
        start= 1,
        end= 2,
        hgnc_symbol="TIMMDC1"
    )

    expect_is(createFullLinkTable(dt, TRUE), "data.table")
    expect_is(createFullLinkTable(dt, FALSE), "data.table")
    expect_is(createFullLinkTable(dt[0]), "data.table")
})
