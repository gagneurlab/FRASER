#'
#' ### MODEL
#'
#' @noRd
fit_autoenc <- function(fds, type=currentType(fds), q_guess=round(ncol(fds)/4),
                        implementation, noiseRatio=0.5, BPPARAM=bpparam(),
                        iterations=3, verbose=FALSE){
    
    message(paste(date(), "; q:", q_guess, "; noise: ", noiseRatio))
    
    # setup object
    currentType(fds) <- type
    verbose_old <- verbose(fds)
    verbose(fds) <- verbose
    on.exit({ verbose(fds) <- verbose_old })
    
    # train AE
    fds <- fit(fds, type=type, q=q_guess, implementation=implementation,
                iterations=iterations, BPPARAM=BPPARAM)
    curLoss <- metadata(fds)[[paste0('loss_', type)]]
    curLoss <- mean(curLoss[,ncol(curLoss)])
    
    return(list(fds=fds, evaluation=curLoss))
}

predict_outliers <- function(fds, type, implementation, BPPARAM){

    fds <- calculatePvalues(fds, type=type, implementation=implementation,
            BPPARAM=BPPARAM)
    fds <- calculatePadjValues(fds, type=type, geneLevel=FALSE, 
            BPPARAM=BPPARAM)

    return(fds)
}

eval_prot <- function(fds, type){
    
    message(date(), ": Calculating AUC-PR ...")
    
    index <- getSiteIndex(fds, type)
    idx   <- !duplicated(index)

    scores <- -as.vector(pVals(fds, type=type)[idx,])

    dt <- cbind(data.table(id=index),
            as.data.table(assay(fds, paste0("trueOutliers_", type))))
    setkey(dt, id)
    labels <- as.vector(vapply(samples(fds), function(i){
        dttmp <- dt[,any(get(i) != 0),by=id]
        setkey(dttmp, id)
        dttmp[J(unique(index)), V1]
    }, FUN.VALUE=logical(length(unique(index))) ) ) + 0

    if(any(is.na(scores))){
        # warning(sum(is.na(scores)), " P-values where NAs.")
        scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
    }
    pr <- pr.curve(scores, weights.class0=labels)
    pr

    return(pr$auc.integral)
}


findEncodingDim <- function(i, fds, type, params, implementation,
                            internalBPPARAM=1, iterations){
    iBPPARAM <- getBPParam(internalBPPARAM)
    
    q_guess    <- params[i, "q"]
    noiseRatio <- params[i, "noise"]
    message(paste(i, ";\t", q_guess, ";\t", noiseRatio))
    implementation <- getHyperOptimCorrectionMethod(implementation)
    
    res_fit <- fit_autoenc(fds=fds, type=type, q_guess=q_guess,
            implementation=implementation, 
            noiseRatio=noiseRatio, BPPARAM=iBPPARAM,
            iterations=iterations)
    res_pvals <- predict_outliers(res_fit$fds, implementation=implementation,
            type=type, BPPARAM=iBPPARAM)
    evals <- eval_prot(res_pvals, type=type)

    return(list(q=q_guess, noiseRatio=noiseRatio, loss=res_fit$evaluation, 
                aroc=evals))
}

#'
#' Find optimal encoding dimension
#'
#' Finds the optimal encoding dimension by either Optimal Hard Thresholding 
#' or injecting artificial splicing outlier ratios while maximizing the 
#' precision-recall curve.
#' 
#' @inheritParams countRNA
#' @inheritParams FRASER
#' @inherit fit
#' @param useOHT If \code{TRUE} (default), Optimal Hard Thresholding (OHT) is 
#' used to estimate the optimal encoding dimension. OHT is only supported for 
#' the Intron Jaccard Index. 
#' @param q_param Vector specifying which values of q should be tested
#' @param noise_param Vector specifying which noise levels should be tested.
#' @param setSubset The size of the subset of the most variable introns that 
#' should be used for the hyperparameter optimization.
#' @param internalThreads The number of threads used internally.
#' @param injectFreq The frequency with which outliers are injected into the 
#' data.
#' @param plot If \code{TRUE}, a plot of the area under the curve and the 
#' model loss for each evaluated parameter combination will be displayed after 
#' the hyperparameter optimization finishes.
#' @param delayed If FALSE, count matrices will be loaded into memory (faster 
#' calculations), otherwise the function works on the delayedMatrix 
#' representations (more memory efficient). The default value depends on the 
#' number of samples in the fds-object.
#' @param ... Additional parameters passed to \code{injectOutliers}.
#' 
#' @return FraserDataSet
#'
#' @examples
#'   # generate data
#'   fds <- makeSimulatedFraserDataSet(m=15, j=20)
#'   fds <- calculatePSIValues(fds)
#'   
#'   # run OHT
#'   fds <- estimateBestQ(fds, type="jaccard", useOHT=TRUE)
#'   
#'   # run hyperparameter optimization
#'   fds <- estimateBestQ(fds, type="jaccard", useOHT=FALSE, q_param=c(2, 5))
#'   
#'   # get estimated optimal dimension of the latent space
#'   bestQ(fds, type="jaccard")
#'   hyperParams(fds, type="jaccard")
#'   
#' @export
estimateBestQ <- function(fds, type="jaccard", useOHT=TRUE, implementation="PCA",
                    q_param=getEncDimRange(fds),
                    noise_param=0, minDeltaPsi=0.1,
                    iterations=5, setSubset=50000, injectFreq=1e-2,
                    BPPARAM=bpparam(), internalThreads=1, plot=TRUE, 
                    delayed=ifelse(ncol(fds) <= 300, FALSE, TRUE), ...){
    
    # Optimal Hard Thresholding (default)
    if (isTRUE(useOHT)){
      if (type=="jaccard"){
        
        # Calculate centered autoencoder input
        x <- as.matrix(x(fds, type=type, all=TRUE, noiseAlpha = NULL, center = TRUE))
        
        # Check for infinite values
        if (any(is.infinite(x))) {
          stop("OHT input matrix contains infinite values.")
        }
        
        message(date(), ": Estimating the optimal encoding dimension using OHT...")
        
        # Transpose x if aspect ratio larger than 1
        if (ncol(x)/nrow(x) > 1){
          x <- t(x)
        } 
        
        # Perform Singular Value Decomposition (SVD) on the matrix 
        # and extract singular values
        sv <- svd(x)$d
        
        # Aspect ratio of the (transposed) matrix
        numRows <- nrow(x)
        numCols <- ncol(x)
        beta <- numCols / numRows
        
        # Compute the optimal w(beta)
        coef <- (optimalSVHTCoef(beta) / sqrt(medianMarchenkoPastur(numCols, numRows)))
        
        # compute cutoff
        cutoff <- coef * median(sv)
        
        # compute and return rank
        if (any(sv > cutoff)){
          latentDim <- max(which(sv > cutoff))
        } 
        else {
          warning(paste("Optimal latent space dimension is smaller than 2. Check your count matrices and",
                        "verify that all samples have the expected number of counts",
                        "For now, the latent space dimension is set to 2.", collapse = "\n"))
          latentDim <- 2
        } 
        cat("Optimal encoding dimension:", latentDim, "\n")
        
        data <- data.table(
          q=latentDim,
          oht=cutoff
          )
      
        hyperParams(fds, type=type) <- data
        #TODO: integrate singular value plot into plotEncDimSearch function
        #if(isTRUE(plot)){
        #  print(plotEncDimSearch(fds, type=type, plotType="sv"))
        #}
        
        return(fds)
      }
      
      else{
        stop(paste("OHT is only supported for the Intron Jaccard Index.",
             "Set useOHT=FALSE to perform a hyperparameter optimization",
             "for other metrics.", collapse = "\n"))
      }
    } 
  
    # Hyperparamter optimization
    else{
      if(isFALSE(needsHyperOpt(implementation))){
          message(date(), ": For correction '", implementation, "' no hyper ",
                  "parameter optimization is needed.")
          data <- data.table(q=NA, noise=0, eval=1, aroc=1)
          hyperParams(fds, type=type) <- data
          return(fds)
      }
  
      #
      # put the most important stuff into memory
      #
      currentType(fds) <- type
      if(isFALSE(delayed)){
          counts(fds, type=type, side="other", HDF5=FALSE) <-
                  as.matrix(counts(fds, type=type, side="other"))
          counts(fds, type=type, side="ofInterest", HDF5=FALSE) <-
                  as.matrix(counts(fds, type=type, side="ofInterest"))
      }
  
      #
      # remove non variable and low abundance junctions
      #
      j2keepVa <- variableJunctions(fds, type, minDeltaPsi)
      j2keepDP <- rowQuantiles(K(fds, type), probs=0.75, drop=FALSE)[,1] >= 10
      j2keep <- j2keepDP & j2keepVa
      message("dPsi filter:", pasteTable(j2keep))
      # TODO fds <- fds[j2keep,,by=type]
  
      # ensure that subset size is not larger that number of introns/splice sites
      setSubset <- pmin(setSubset, nrow(K(fds, type)))
      
      optData <- data.table()
      for(nsub in setSubset){
          # subset for finding encoding dimensions
          # most variable functions + random subset for decoder
          exMask <- subsetKMostVariableJunctions(fds[j2keep,,by=type], type, nsub)
          j2keep[j2keep==TRUE] <- exMask
  
          # keep n most variable junctions + random subset
          j2keep <- j2keep | sample(c(TRUE, FALSE), length(j2keep), replace=TRUE,
                  prob=c(nsub/length(j2keep), 1 - nsub/length(j2keep)))
          message("Exclusion matrix: ", pasteTable(j2keep))
  
          # make copy for testing
          fds_copy <- fds
          dontWriteHDF5(fds_copy) <- TRUE
          featureExclusionMask(fds_copy) <- j2keep
          fds_copy <- fds_copy[j2keep,,by=type]
          currentType(fds_copy) <- type
  
          # inject outliers
          fds_copy <- injectOutliers(fds_copy, type=type, freq=injectFreq, ...)
          
          if(sum(getAssayMatrix(fds_copy, type=type, "trueOutliers") != 0) == 0){
              warning(paste0("No outliers could be injected so the ", 
                              "hyperparameter optimization could not run. ", 
                              "Possible reason: too few junctions in the data."))
              return(fds)
          }
  
          # remove unneeded blocks to save memory
          a2rm <- paste(sep="_", c("originalCounts", "originalOtherCounts"),
                  rep(psiTypes, 2))
          for(a in a2rm){
              assay(fds_copy, a) <- NULL
          }
          metadata(fds_copy) <- list()
          gc()
  
          # reset lost important values
          currentType(fds_copy) <- type
          dontWriteHDF5(fds_copy) <- TRUE
  
          # run hyper parameter optimization
          params <- expand.grid(q=q_param, noise=noise_param)
          message(date(), ": Run hyper optimization with ", nrow(params), 
                  " options.")
          res <- bplapply(seq_len(nrow(params)), findEncodingDim, fds=fds_copy, 
                          type=type,
                          iterations=iterations, params=params, 
                          implementation=implementation,
                          BPPARAM=BPPARAM, 
                          internalBPPARAM=internalThreads)
  
          data <- data.table(
              q=vapply(res, "[[", "q", FUN.VALUE=numeric(1)),
              noise=vapply(res, "[[", "noiseRatio", FUN.VALUE=numeric(1)),
              nsubset=nsub,
              eval=vapply(res, "[[", "loss", FUN.VALUE=numeric(1)),
              aroc=vapply(res, "[[", "aroc", FUN.VALUE=numeric(1)))
  
          optData <- rbind(optData, data)
      }
  
      hyperParams(fds, type=type) <- optData
      if(isTRUE(plot)){
          print(plotEncDimSearch(fds, type=type, plotType="auc"))
      }
      return(fds)
    }
}

#'
#' Calculate the OHT coefficient 
#' 
#' @noRd
optimalSVHTCoef <- function(beta){ 
  # Calculate lambda(beta)
  sqrt(2 * (beta + 1) + (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1)))
}

#'
#' Calculate the median of the Marchenko-Pastur distribution
#' 
#' Formulas are derived from a robust estimator for the noise 
#' parameter gamma in the model Z~ = Z + gamma * E.
#' Gamma(Z~) = sigma / sqrt(n * mu)) with mu: median of the Marcenko-Pastur distribution.
#' More detailed explanations can be found in Gavish and Donoho 2014.
#' 
#' @noRd
medianMarchenkoPastur <- function(ncol, nrow){ 
  # Compute median of Marchenko-Pastur distribution
  beta <- ncol / nrow
  betaMinus <- (1 - sqrt(beta))^2
  betaPlus <- (1 + sqrt(beta))^2
  
  # Initialize range for upper integral boundary
  lobnd <- copy(betaMinus) 
  hibnd <- copy(betaPlus)
  
  
  while ((hibnd - lobnd) > 0.001){ # Iterate until convergence
    x <- seq(lobnd, hibnd, length.out = 10) # Set range of values for upper integral boundary
    y <- rep(0, length(x))
    for (numCols in 1:length(x)){
      # Approximate integral using Gauss-Kronrod Quadrature
      y[numCols] <- pracma::quadgk(RMTstat::dmp, a=betaMinus, b=x[numCols], ndf=nrow, pdim=ncol)
    }  
    
    # Set new boundaries for x that yield the closest results to 0.5
    if (any(y < 0.5)){
      lobnd = max(x[y < 0.5])
    }
    
    if (any(y > 0.5)){
      hibnd = min(x[y > 0.5])
    }
  }
  # If hibnd and lobnd are similar enough, return their mean
  return((hibnd + lobnd) / 2.0)
}

#' Get default range of latent space dimensions to test during hyper param opt
#' @noRd
getEncDimRange <- function(fds, mp=3){
    # Get range for latent space dimension
    a <- 2 
    b <- min(ncol(fds), nrow(fds)) / mp   # N/mp
    
    maxSteps <- 12
    if(mp < 6){
        maxSteps <- 15
    }
    
    Nsteps <- min(maxSteps, b)
    pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
    return(pars_q)
}
