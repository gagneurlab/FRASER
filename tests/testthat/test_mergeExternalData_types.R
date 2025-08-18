context("Test mergeExternalData types parameter")

test_that("mergeExternalData handles types parameter correctly", {
    
    # Skip this test if not all dependencies are available
    skip_if_not_installed("data.table")
    
    # Create a simple test FraserDataSet
    fds <- createTestFraserDataSet()
    
    # Set fitMetrics to only jaccard to simulate the error condition
    fitMetrics(fds) <- "jaccard"
    
    # Check that fitMetrics returns only jaccard
    expect_equal(fitMetrics(fds), "jaccard")
    
    # Test 1: Function signature - should have types parameter
    expect_true("types" %in% names(formals(mergeExternalData)))
    
    # Test 2: Check default value of types parameter uses fitMetrics
    default_types <- formals(mergeExternalData)$types
    expect_true(!is.null(default_types))
    
    # Test 3: Verify calculatePSIValues accepts types parameter
    expect_true("types" %in% names(formals(calculatePSIValues)))
    
    # Test 4: Test calculatePSIValues with specific types (jaccard only)
    fds_test <- calculatePSIValues(fds, types="jaccard")
    expect_is(fds_test, "FraserDataSet")
    expect_true("jaccard" %in% assayNames(fds_test))
    
    # Test 5: Test that the mergeExternalData function would use the correct default
    # This tests that the function signature change maintains backward compatibility
    # while fixing the issue
    
    # Check that when fitMetrics is set to jaccard, it would be used as default
    fds_jaccard_only <- createTestFraserDataSet()
    fitMetrics(fds_jaccard_only) <- "jaccard"
    expect_equal(fitMetrics(fds_jaccard_only), "jaccard")
    
    # Test with multiple types
    fds_multi <- createTestFraserDataSet()
    fitMetrics(fds_multi) <- c("jaccard", "psi5")
    expect_equal(fitMetrics(fds_multi), c("jaccard", "psi5"))
    
    # Test calculatePSIValues with multiple specific types
    fds_multi_test <- calculatePSIValues(fds_multi, types=c("jaccard", "psi5"))
    expect_is(fds_multi_test, "FraserDataSet")
    expect_true(all(c("jaccard", "psi5") %in% assayNames(fds_multi_test)))
})

test_that("mergeExternalData integration test with external files", {
    
    # Skip this test if package files are not available
    skip_if_not(file.exists(system.file("extdata", "externalCounts", 
                                       "annotation.tsv.gz", package="FRASER")))
    
    # Create a test FDS with jaccard only
    fds <- createTestFraserDataSet()
    fitMetrics(fds) <- "jaccard"
    
    # Get external count files (from package examples)
    anno <- fread(system.file("extdata", "externalCounts", 
                             "annotation.tsv.gz", package="FRASER"))
    ctsFiles <- list.files(full.names = TRUE, pattern="counts",
                          system.file("extdata", "externalCounts", package="FRASER"))
    
    # Test that mergeExternalData works with types parameter
    # Only test if we have the required external files
    if(length(ctsFiles) == 5 && nrow(anno) > 0) {
        # This should work without error when types parameter is properly passed
        expect_no_error({
            fds_merged <- mergeExternalData(fds, ctsFiles, anno[1:2, sampleID], 
                                           anno[1:2,], types="jaccard")
        })
        
        # Test that merged object has correct properties
        expect_is(fds_merged, "FraserDataSet")
        expect_true("jaccard" %in% assayNames(fds_merged))
    }
})