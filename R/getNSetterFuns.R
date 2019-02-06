
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
featureExclusionMask <- function(fds, type=currentType(fds)){
    ans <- rep(FALSE, nrow(mcols(fds, type=type)))
    if(paste0('featureExclude_', type) %in% colnames(mcols(fds, type=type))){
        ans <- mcols(fds, type=type)[[paste0('featureExclude_', type)]]
    }
    # TODO names(ans, type=type) <- rownames(fds, type=type)
    return(ans)
}

#' @rdname featureExclusionMask
#' @export "featureExclusionMask<-"
`featureExclusionMask<-` <- function(fds, value, type=currentType(fds)){
    if(isScalarLogical(value)){
        value <- rep(value, ncol(mcols(fds, type=type)))
    }
    mcols(fds, type=type)[[paste0('featureExclude_', type)]] <- value
    return(fds)
}

K <- function(fds, type=currentType(fds)){
  K <- counts(fds, type=type, side="ofInterest")
  return(K);
}

N <- function(fds, type=currentType(fds)){
  N <- K(fds, type=type) + counts(fds, type=type, side="other")
  return(N);
}

x <- function(fds, type=currentType(fds), all=FALSE){
  K <- K(fds, type=type)
  N <- N(fds, type=type)

  # compute logit ratio with pseudocounts
  x <- t((K + 0.5)/(N + 1))
  x <- qlogis(x)

  if(isFALSE(all)){
      x = x[,featureExclusionMask(fds, type=type)]
  }
  return(x)
}

H <- function(fds, type=currentType(fds)){
  x(fds, all=FALSE, type=type) %*% E(fds, type=type)
}

`D<-` <- function(fds, value, type=currentType(fds)){
  if(!is.matrix(value)){
    value <- matrix(value, nrow=nrow(fds))
  }
  metadata(fds)[[paste0('D_', type)]] <- value
  return(fds)
}

D <- function(fds, type=currentType(fds)){
  return(metadata(fds)[[paste0('D_', type)]])
}

`E<-` <- function(fds, value, type=currentType(fds)){
  if(!is.matrix(value)){
    value <- matrix(value, nrow=sum(featureExclusionMask(fds, type=type)))
  }
  metadata(fds)[[paste0('E_', type)]] <- value
  return(fds)
}

E <- function(fds, type=currentType(fds)){
  return(metadata(fds)[[paste0('E_', type)]])
}

`b<-` <- function(fds, value, type=currentType(fds)){
  mcols(fds, type=type)[[paste0('b_', type)]] <- value
  return(fds)
}

b <- function(fds, type=currentType(fds)){
  return(mcols(fds, type=type)[[paste0('b_', type)]])
}

`rho<-` <- function(fds, value, type=currentType(fds)){
  mcols(fds, type=type)[[paste0('rho_', type)]] <- value
  return(fds)
}

rho <- function(fds, type=currentType(fds)){
  return(mcols(fds, type=type)[[paste0('rho_', type)]])
}

predictMu <- function(fds, type=currentType(fds)){
  y <- predictY(fds, type=type)
  mu <- predictMuCpp(y)
  return(t(mu))
}

predictY <- function(fds, type=currentType(fds)){
  D <- D(fds, type=type)
  b <- b(fds, type=type)
  H <- H(fds, type=type)

  y <- predictYCpp(as.matrix(H), D, b)

  return(t(y))
}

zScores <- function(fds, type=currentType(fds)){
  return(assay(fds, paste0('zScores_', type=type)))
}

`zScores<-` <- function(fds, value, type=currentType(fds), ...){
  if(!is.matrix(value)){
    value <- matrix(value, ncol=ncol(fds))
  }
  assay(fds, paste0('zScores_', type), ...) <- value
  return(fds)
}

pVals <- function(fds, type=currentType(fds)){
  return(assay(fds, paste0('pvalues_', type)))
}

`pVals<-` <- function(fds, value, type=currentType(fds), ...){
  if(!is.matrix(value)){
    value <- matrix(value, ncol=ncol(fds))
  }
  assay(fds, paste0('pvalues_', type), ...) <- value
  return(fds)
}

predictedMeans <- function(fds, type=currentType(fds)){
    return(assay(fds, paste0('predictedMeans_', type)))
}

`predictedMeans<-` <- function(fds, value, type=currentType(fds), ...){
    if(!is.matrix(value)){
        value <- matrix(value, ncol=ncol(fds))
    }
    assay(fds, paste0('predictedMeans_', type), ...) <- value
    return(fds)
}

currentType <- function(fds){
    return(metadata(fds)[['currentType']])
}

`currentType<-` <- function(fds, value){
    stopifnot(isScalarCharacter(whichPSIType(value)))
    metadata(fds)[['currentType']] <- whichPSIType(value)
    return(fds)
}
