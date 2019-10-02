context("Testing the FraseR shiny app")

# testSkipShiny <- FALSE
# 
# tmpTry <- try(silent=TRUE, {
#     remDr <- remoteDriver()
#     remDr$open(silent = TRUE)
#     appURL <- "http://ouga03:3843/"
# })
#  
# if(class(tmpTry) == "try-error"){
#     message("Skip shiny app test since we can not connect to server.")
#     testSkipShiny <- TRUE
# }
# 
# test_that("can connect to app", {
#     if(testSkipShiny) return(TRUE)
# 
#     remDr$navigate(appURL)
#     appTitle <- remDr$getTitle()[[1]]
#     expect_equal(appTitle, "FraseR")
# })
# 
# try({remDr$close()}, silent=TRUE)
# 
