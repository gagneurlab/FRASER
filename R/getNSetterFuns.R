
K <- function(fds){
  K <- counts(fds, type="psi3", side="ofInterest")
  return(K);
}

N <- function(fds){
  N <- K(fds) + counts(fds, type="psi3", side="other")
  return(N);
}

x <- function(fds){
  K <- K(fds)
  N <- N(fds)
  
  # compute logit ratio with pseudocounts
  x <- t((K + 0.5)/(N + 1))
  x <- qlogis(x)
  
  return(as.matrix(x))
}

H <- function(fds){
  x(fds) %*% E(fds)
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
    value <- matrix(value, nrow=nrow(fds))
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

predictMu <- function(fds){
  D <- D(fds)
  b <- b(fds)
  H <- H(fds)
  
  mu <- predictMuCpp(H, D, b)
  
  return(t(mu))
}