#'
#' Update E step for the autoencoder fit
#' 
#' @noRd
updateE <- function(fds, control, BPPARAM, verbose){
  par <- as.vector(E(fds))
  D <- D(fds)
  k <- as.matrix(counts(fds, type="psi3", side="ofInterest"))
  o <- as.matrix(counts(fds, type="psi3", side="other"))
  n <- k + o
  x <- x(fds)
  b <- b(fds)
  rho <- rho(fds)
  
  control[['maxit']] <- 50
  #control[['parscale']] <- rep(4, length(par))
  
  fit <- optim(par, fn=truncNLL_e, gr=truncGrad_e, 
               x=x, D=D, k=k, n=n, rho=rho, b=b, 
               method="L-BFGS-B", control=control)
  
  # Check that fit converged
  if(isTRUE(verbose) & fit$convergence != 0){
    print(paste('Update E did not converge: ', fit$message))
  }
  
  # update ods
  E(fds) <- fit$par
  
  return(fds)
}
