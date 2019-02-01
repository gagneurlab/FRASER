
#'
#' Feature exclusion
#'
#' To remove certain junctions from being used in the train step of the encoding dimension
#' we can set the \code{featureExclusion} vector to \code{FALSE}. This can be helpfull if
#' we have local linkage between features which we do not want to model by the autoencoder.
#'
#' @param ods An OutriderDataSet object
#' @param value A logical vector of the length of the features. If \code{TRUE},
#'             the corresponding feature will be excluded from the encoding dimension
#'             fit.
#' @return The exclusion vector
#'
#' @name featureExclusionMask
#' @rdname featureExclusionMask
#' @aliases featureExclusionMask, `featureExclusionMask<-`
#'
#' @examples
#' ods <- makeExampleFraseRDataSet()
#' featureExclusionMask(fds) <- sample(c(FALSE, TRUE), nrow(fds), replace=TRUE)
#'
#' featureExclusionMask(fds)
#'
#' @export featureExclusionMask
#' @export "featureExclusionMask<-"
featureExclusionMask <- function(fds){
    ans <- rep(FALSE, ncol(fds))
    if('featureExclude' %in% colnames(mcols(fds))){
        ans <- mcols(fds)[['featureExclude']]
    }
    names(ans) <- rownames(fds)
    return(ans)
}

#' @rdname featureExclusionMask
#' @export "featureExclusionMask<-"
`featureExclusionMask<-` <- function(fds, value){
    if(isScalarLogical(value)){
        value <- rep(value, ncol(fds))
    }
    mcols(fds)[['featureExclude']] <- value
    return(fds)
}

K <- function(fds){
  K <- counts(fds, type="psi3", side="ofInterest")
  return(K);
}

N <- function(fds){
  N <- K(fds) + counts(fds, type="psi3", side="other")
  return(N);
}

x <- function(fds, all=FALSE){
  K <- K(fds)
  N <- N(fds)

  # compute logit ratio with pseudocounts
  x <- t((K + 0.5)/(N + 1))
  x <- qlogis(x)

  if(isFALSE(all)){
      x = x[,featureExclusionMask(fds)]
  }
  return(as.matrix(x))
}

H <- function(fds, all=FALSE){
  x(fds, all) %*% E(fds)
}

`D<-` <- function(fds, value){
  if(!is.matrix(value)){
    value <- matrix(value, nrow=nrow(fds))
  }
  metadata(fds)[['D']] <- value
  return(fds)
}

D <- function(fds){
  return(metadata(fds)[['D']])
}

`E<-` <- function(fds, value){
  if(!is.matrix(value)){
    value <- matrix(value, nrow=sum(featureExclusionMask(fds)))
  }
  metadata(fds)[['E']] <- value
  return(fds)
}

E <- function(fds){
  return(metadata(fds)[['E']])
}

`b<-` <- function(fds, value){
  mcols(fds)[['b']] <- value
  return(fds)
}

b <- function(fds){
  return(mcols(fds)[['b']])
}

`rho<-` <- function(fds, value){
  mcols(fds)[['rho']] <- value
  return(fds)
}

rho <- function(fds){
  return(mcols(fds)[['rho']])
}

predictMu <- function(fds, all=FALSE){
  y <- predictY(fds, all)
  mu <- predictMuCpp(y)
  return(mu)
}

predictY <- function(fds, all=FALSE){
  D <- D(fds)
  b <- b(fds)
  H <- H(fds, all)

  y <- predictYCpp(H, D, b)

  return(t(y))
}

zScores <- function(fds){
  return( metadata(fds)[['zScores']] )
}

`zScores<-` <- function(fds, value){
  if(!is.matrix(value)){
    value <- matrix(value, nrow=nrow(fds))
  }
  metadata(fds)[['zScores']] <- value
  return(fds)
}

pVals <- function(fds){
  return( metadata(fds)[['pvalues']] )
}

`pVals<-` <- function(fds, value){
  if(!is.matrix(value)){
    value <- matrix(value, nrow=nrow(fds))
  }
  metadata(fds)[['pvalues']] <- value
  return(fds)
}
