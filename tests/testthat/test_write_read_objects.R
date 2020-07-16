context("testing HDF5 write/read functionality")

test_that("save/load/move fds object", {
    # temp folder
    dir1 <- tempfile("fraser_test_1_")
    dir2 <- tempfile("fraser_test_2_")
    on.exit(unlink(dir1, force=TRUE, recursive=TRUE))
    on.exit(unlink(dir2, force=TRUE, recursive=TRUE))
    
    # create test object (force to write hdf5)
    options(FRASER.maxJunctionsNoHDF5=10)
    on.exit(options(FRASER.maxJunctionsNoHDF5=1000))
    fds <- makeSimulatedFraserDataSet(workingDir=dir1, j=11)
    fds <- saveFraserDataSet(fds)
    
    # load from created folder
    fds_l1 <- loadFraserDataSet(file=file.path(
            dir1, "savedObjects", "Data_Analysis", "fds-object.RDS"))
    expect_equivalent(fds_l1, fds)
    expect_match(path(assay(fds, "rawCountsJ")), dir1)
    
    # rename and check if still the same
    renameFile(dir1, dir2)
    fds_l2 <- loadFraserDataSet(file=file.path(
            dir2, "savedObjects", "Data_Analysis", "fds-object.RDS"))
    expect_equivalent(fds_l2, fds)
    expect_match(path(assay(fds_l2, "rawCountsJ")), dir2)
})

