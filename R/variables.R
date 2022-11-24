#'
#' Available splice metrics
#'
#' @examples 
#'   # to show available splice metrics:
#'   psiTypes_avail
#' 
#'   # to show splice metrics selected to be fitted:
#'   psiTypes
#' 
#' @rdname psiTypes
#' @export
psiTypes_avail <- c("jaccard", "psi5", "psi3", "theta")
names(psiTypes_avail) <- c("Intron Jaccard Index", "psi5", "psi3", "theta")

#' @describeIn psiTypes Splice metrics that are run by default
#' @export
psiTypes <- c("jaccard")
names(psiTypes) <- c("Intron Jaccard Index")
