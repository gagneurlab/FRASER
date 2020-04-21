#
# helper scripts for the testing step
#
getFraser <- function(clean=FALSE){
    fds <- createTestFraserDataSet(rerun=clean)
    fds
}
