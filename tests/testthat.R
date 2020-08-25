library(testthat)
library(FRASER)

# to speed up the testing on windows do it in serial mode
if(.Platform$OS.type != "unix") {
    register(SerialParam())
}

test_check("FRASER")
