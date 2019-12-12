#'
#' Create example data sets for FraseR
#'
#' Creates an example data set from a file or simulates a data set based
#' on  random split read counts following a beta-binomial distribution with
#' injected outliers.
#' 
#' @return FraseRDataSet
#'
#' @examples
#'   fds <- makeExampleFraseRDataSet()
#'   fds
#'
#' @export
makeExampleFraseRDataSet <- function(){
    # TODO remove create and change it to makeExampleFraseRDataSet
    createTestFraseRDataSet()
}
