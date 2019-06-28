#'
#' Update D function
#'
#' @noRd
updateD <- function(fds, type, lambda, control, BPPARAM, verbose, nrDecoderBatches=5, weighted=FALSE){
  D <- D(fds)
  b <- b(fds)
  H <- H(fds, noiseAlpha=currentNoiseAlpha(fds))
  k <- K(fds)
  n <- N(fds)
  rho <- rho(fds)

  # run fits
  fitls <- bplapply(seq_len(nrow(k)), singleDFit, D=D, b=b, k=k, n=n, H=H,
                    rho=rho, lambda=lambda, control=control, weighted=weighted,
                    nSamples=ncol(fds), nrBatches=nrDecoderBatches, BPPARAM=BPPARAM)

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

singleDFit <- function(i, D, b, k, n, H, rho, lambda, control, nSamples, nrBatches, weighted, ...){
  pari <- c(b[i], D[i,])
  ki <- k[i,]
  ni <- n[i,]
  rhoi <- rho[i]

  if(nrBatches == 1){

    fit <- optim(pari, fn=truncWeightedNLL_db, gr=truncWeightedGrad_db,
                 H=H, k=ki, n=ni, rho=rhoi, lambda=lambda, control=control,
                 weighted=weighted, method="L-BFGS-B")

  }
  else { # cross-validation

    subsets <- split(sample(1:nSamples), rep(1:nrBatches, length = nSamples))
    fitls <- sapply(subsets, function(subset){
      ksub <- ki[-subset]
      nsub <- ni[-subset]
      Hsub <- H[-subset,]

      fit <- optim(pari, fn=truncWeightedNLL_db, gr=truncWeightedGrad_db,
                   H=Hsub, k=ksub, n=nsub, rho=rhoi, lambda=lambda, control=control,
                   weighted=weighted, method="L-BFGS-B")
      return(fit)
    })

    pars <- colMedians(do.call(rbind, fitls['par',]))

    fit <- fitls[,1]
    fit$par         <- pars
    fit$value       <- mean(unlist(fitls['value',]))
    fit$counts      <- round(colMedians(do.call(rbind, fitls['counts',])))
    fit$convergence <- sum(unlist(fitls['convergence',]))

  }

  return(fit)
}

