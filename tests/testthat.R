library(testthat)
library(FRASER)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# to speed up the testing on windows do it in serial mode
if(.Platform$OS.type != "unix") {
    register(SerialParam())
}

set.seed(42)
test_check("FRASER")
