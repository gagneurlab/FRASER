#
# helper scripts for the testing step
#
getFraseR <- function(clean=FALSE){
    if(!"initBPPARAM" %in% ls()){
        if(.Platform$OS.type != "windows"){
            register(MulticoreParam(min(4, bpworkers())))
        } else {
            register(SerialParam())
        }
    }
    initBPPARAM <- TRUE
    assign("initBPPARAM", initBPPARAM, envir=parent.env(environment()))
    
    fds <- createTestFraseRDataSet(rerun=clean)
    fds
}
