context("Test merge fds function")

test_that("mergeFDS", {
  fds <- createTestFraserDataSet()
  
  fds_1 <- fds[, 1]
  fds_2 <- fds[, c(2,3)]
  
  merged_fds <- mergeFDS(fds_1, fds_2)
  
  # Check number of samples
  expect_equal(ncol(merged_fds), ncol(fds))
  
  # Check number of junctions (features)
  expect_equal(nrow(merged_fds), nrow(fds))
  
  # Check number of splice sites
  expect_equal(nrow(assay(merged_fds, "rawCountsSS")), nrow(assay(fds, "rawCountsSS")))  
})