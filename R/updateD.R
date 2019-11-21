#'
#' Update D function
#'
#' @noRd
updateD <- function(fds, type, lambda, control, BPPARAM, verbose,
                  nrDecoderBatches=5, weighted=FALSE, fraction=0.75,
                  multiRho=FALSE){
  D <- D(fds)
  b <- b(fds)
  H <- H(fds, noiseAlpha=currentNoiseAlpha(fds))
  k <- K(fds)
  n <- N(fds)
  rho <- rho(fds)

  weights <- matrix(1, nrow=nrow(k), ncol=ncol(k))
  if(isTRUE(weighted)){
    weights <- calcFraseRWeights(fds, type)
    weights(fds, type, HDF5=FALSE) <- weights
    message(sum(c(weights) != 1))
  }

  # run fits
  fitls <- bplapply(seq_len(nrow(k)), singleDFit, D=D, b=b, k=k, n=n, H=H,
                    rho=rho, lambda=lambda, multiRho=multiRho, control=control,
                    weights=weights, nSamples=ncol(fds), fraction=fraction,
                    nrBatches=nrDecoderBatches, BPPARAM=BPPARAM)

  # extract infos
  parMat <- vapply(fitls, '[[', double(ncol(D) + 1), 'par')
  mcols(fds, type=type)[paste0('FitDMessage_', type)] <-
      vapply(fitls, '[[', 'text', 'message')
  mcols(fds, type=type)[,paste0('NumConvergedD_', type)] <-
    mcols(fds, type=type)[,paste0('NumConvergedD_', type)] + grepl(
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

singleDFit <- function(i, D, b, k, n, H, rho, lambda, control, fraction,
                  nrBatches, weights, nSamples, multiRho, ...){
  pari <- c(b[i], D[i,])
  ki <- k[i,]
  ni <- n[i,]
  rhoi <- rho[i]
  wi <- weights[i,]

  if(multiRho){
    rhoi <- c(rhoi/2, rhoi, (1+rhoi)/2)
  }

  # use full set if only one batch is used
  if(nrBatches == 1){
    fraction <- 1
  }

  fitls <- lapply(seq_len(nrBatches), function(subset){
      subset <- sample(seq_len(nrow(H)), ceiling(nrow(H)*fraction))
      ksub <- ki[subset]
      nsub <- ni[subset]
      Hsub <- H[subset,]
      wsub <- wi[subset]

      fit <- lapply(rhoi, function(r){
          optim(pari, fn=truncWeightedNLL_db, gr=truncWeightedGrad_db,
                   H=Hsub, k=ksub, n=nsub, rho=r, lambda=lambda, w=wsub,
                   method="L-BFGS-B", control=control)
      })
      fit[[which.min(vapply(fit, "[[", "value", FUN.VALUE=numeric(1)))]]
  })

  # pars <- colMedians(do.call(rbind, fitls['par',]))
  pars <- rowMedians(vapply(fitls, "[[", "par", FUN.VALUE=pari))

  fit <- fitls[[1]]
  fit$par         <- pars
  fit$value       <- mean(vapply(fitls, "[[", 'value', FUN.VALUE=numeric(1)))
  fit$counts      <- round(rowMedians(vapply(fitls, "[[", 'counts', 
                                             FUN.VALUE=numeric(2))))
  fit$convergence <- sum(vapply(fitls, "[[", 'convergence', 
                                FUN.VALUE=numeric(1)))

  return(fit)
}

