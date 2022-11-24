#'
#' Available splice metrics
#'
#' @examples 
#'   # to show available splice metrics:
#'   psiTypes_avail
#' 
#' @export
psiTypes_avail <- c("psi5", "psi3", "theta", "jaccard")
names(psiTypes_avail) <- c("psi5", "psi3", "theta", "Intron Jaccard Index")

#'
#' Splice metrics that are run by default
#'
#' @examples 
#'   # to show splice metrics selected to be fitted:
#'   psiTypes
#' 
#' @export
psiTypes <- c("jaccard")
names(psiTypes) <- c("Intron Jaccard Index")

