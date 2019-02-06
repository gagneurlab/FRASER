#'
#' Update E step for the autoencoder fit
#'
#' @noRd
updateE <- function(fds, control, BPPARAM, verbose){
  par <- as.vector(E(fds))
  D <- D(fds)
  k <- K(fds)
  n <- N(fds)
  x <- x(fds)
  b <- b(fds)
  rho <- rho(fds)

  control[['maxit']] <- 50

  fit <- optim(par, fn=truncNLL_e, gr=truncGrad_e,
               x=as.matrix(x), D=D, k=as.matrix(k), n=as.matrix(n), rho=rho, b=b,
               method="L-BFGS-B", control=control)

  # Check that fit converged
  if(isTRUE(verbose) & fit$convergence != 0){
    print(paste('Update E did not converge: ', fit$message))
  }

  # update ods
  E(fds) <- fit$par

  return(fds)
}
