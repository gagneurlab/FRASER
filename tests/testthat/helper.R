#
# helper scripts for the testing step
#
getFraseR <- function(clean=FALSE){
    fds <- createTestFraseRDataSet(rerun=clean)
    fds
}
