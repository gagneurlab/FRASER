#'
#' Update D function
#'
#' @noRd
updateD <- function(fds, type, lambda, control, BPPARAM, verbose){
  D <- D(fds)
  b <- b(fds)
  H <- H(fds, noiseAlpha=currentNoiseAlpha(fds))
  k <- K(fds)
  n <- N(fds)
  rho <- rho(fds)

  # run fits
  fitls <- bplapply(seq_len(nrow(k)), singleDFit, D=D, b=b, k=k, n=n, H=H,
                    rho=rho, lambda=lambda, control=control, BPPARAM=BPPARAM)

  # extract infos
  parMat <- vapply(fitls, '[[', double(ncol(D) + 1), 'par')
  mcols(fds, type=type)[paste0('FitDMessage_', type)] <- vapply(fitls, '[[', 'text', 'message')
  mcols(fds, type=type)[,paste0('NumConvergedD_', type)] <- mcols(fds, type=type)[,paste0('NumConvergedD_', type)] + grepl(
    "CONVERGENCE: REL_REDUCTION_OF_F .. FACTR.EPSMCH",
    mcols(fds, type=type)[,paste0('FitDMessage_', type)])

  if(isTRUE(verbose)){
    print(table(mcols(fds, type=type)[,paste0('FitDMessage_', type)]))
  }

  # update b and D
  b(fds) <- parMat[1,]
  D(fds) <- t(parMat)[,-1]
  metadata(fds)[[paste0('Dfits_', type)]] <- fitls

  return(fds)
}


singleDFit <- function(i, D, b, k, n, H, rho, lambda, control, ...){
  pari <- c(b[i], D[i,])
  ki <- k[i,]
  ni <- n[i,]
  rhoi <- rho[i]

  fit <- optim(pari, fn=truncNLL_db, gr=truncGrad_db,
               H=H, k=ki, n=ni, rho=rhoi, lambda=lambda, control=control,
               method="L-BFGS-B")

  return(fit)
}
