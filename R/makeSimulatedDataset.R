
#' 
#' Create an simulated example data set for FraseR
#' 
#' Simulates a data set based on random counts following a 
#' negative binomial distribution with injected outliers with 
#' a fixed z score away from the mean of the gene.
#' 
#' @param j Number of simulated junctions 
#' @param m Number of simulated samples
#' @param freq Frequency of in-silico outliers
#' @param zScore Absolute z score of in-silico outliers (default 6).
#' @param inj Determines whether counts are injected with the strategy 
#'            ('both', 'low', 'high'), default is 'both'.
#' @param q number of simulated latent variables.
#' 
#' @return An FraserDataSet containing an example dataset based on 
#'         simulated data
#' 
#' @examples
#' # A generic dataset 
#' fds1 <- makeSimulatedFraserDataSet()
#' fds1
#' 
#' # A generic dataset with specificed sample size and injection method
#' fds2 <- makeSimulatedFraserDataSet(m=100, j=500, inj='low')
#' fds2
#' 
#' @export
makeSimulatedFraserDataSet <- function(m=200, j=1000, q=10, freq=1E-3, zScore=4,
                                       inj=c('both', 'low', 'high')){
  
  # nMean <- ceiling(rlnorm(j, 5)) + 10     
  # theta <- rlnorm(j, 4)                         
  nMean <- rnorm(j, 348, 10)                     # nbinom means for n (mean(N) on kremer dataset = 347.62, median=37)   
  theta <- rnorm(j, 0.27, 0.01)                  # nbinom dispersion (~0.27 on N from kremer dataset)
  
  # simulate counts at each junction
  n <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
  nonSplit <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
  
  
  #
  # Create FraseR data set
  #
  sampleIDs <- paste0("sample", seq_len(m))
  anno <- data.table(sampleID = sampleIDs, bamFile=rep(NA, m))
  fds <- FraseRDataSet(colData=anno)
  
  # put in n as rawcountsJ first so it doesn't complain later when assinging k to it
  junctionData <- SummarizedExperiment(
    colData=colData(fds),
    assays=list(rawCountsJ=n),
    rowRanges=GRanges(seqnames=rep("chr1", j), ranges=IRanges(start=rep(0, j), end= rep(1, j)))
  )
  nonSplitData <- SummarizedExperiment(
    colData=colData(fds),
    assays=list(rawCountsSS=nonSplit),
    rowRanges=GRanges(seqnames=rep("chr1", j), ranges=IRanges(start=rep(0, j), end= rep(1, j)))
  )
  fds <- new("FraseRDataSet", 
             junctionData,
             name            = name(fds),
             method          = method(fds),
             parallel        = parallel(fds),
             bamParam        = scanBamParam(fds),
             strandSpecific  = strandSpecific(fds),
             workingDir      = workingDir(fds),
             nonSplicedReads = nonSplitData
  )
  metadata(fds)[['optimalEncDim']]  <- q
  metadata(fds)[['encDimTable']]    <- data.table(
    encodingDimension=q, evaluationLoss=1, evalMethod='simulation')

  
  #
  # Simulate count matrices for all psi types
  #
  for(type in c("psi3", "psiSite")){
    
    rho <- abs(rnorm(j, mean=0.0001, sd=0.05))     # betaBin dispersion
    sdVec <- rep(0.5, m)                           # sd for H matrix
  
    #
    # Simulate covariates.
    #
    H_true <- matrix(rnorm(m*q), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, sd=sdVec), nrow=j, ncol=q)
    y_true <- D_true %*% t(cbind(H_true))
    mu     <- predictMuCpp(y_true)
  
    # Simulate counts
    if(type== "psi3"){
      k <-matrix(rbetabinom(j*m, size=n, prob=mu, rho=rho), nrow=j, ncol=m)
    } else{
      k <-matrix(rbetabinom(j*m, size=nonSplit, prob=mu, rho=rho), nrow=j, ncol=m)
    }
    mode(n) <- 'integer'
    mode(k) <- 'integer'
  
    psi <- (k + pseudocount())/(n + 2*pseudocount())
    l2fc <- log2fc(psi, mu)
    datasd <- rowSds(l2fc)
    datasd <- matrix(datasd, nrow=j, ncol=m)
    lmu <- log2(mu)
    
    setAssayMatrix(fds=fds, name="truePSI", type=type) <- mu 
    setAssayMatrix(fds=fds, name="trueSd", type=type)  <- datasd
    
    # needed so that subsetting the fds works later
    mcols(fds, type=type)[["startID"]] <- 1:nrow(mcols(fds, type=type))
    mcols(fds, type=type)[["endID"]] <- 1:nrow(mcols(fds, type=type))
    mcols(fds, type=type)[,"trueRho"] <- rho
    
    # for now: same values for psi3 and psi5
    if(type== "psi3"){
      setAssayMatrix(fds=fds, name="truePSI", type="psi5") <- mu 
      setAssayMatrix(fds=fds, name="trueSd", type="psi5")  <- datasd
      
      # needed so that subsetting the fds works later
      mcols(fds, type="psi5")[["startID"]] <- 1:nrow(mcols(fds, type=type))
      mcols(fds, type="psi5")[["endID"]] <- 1:nrow(mcols(fds, type=type))
      mcols(fds, type="psi5")[,"trueRho"] <- rho
    } else{
      # needed so that subsetting the fds works later
      mcols(fds, type=type)[["spliceSiteID"]] <- 1:nrow(mcols(fds, type=type))
    }
    
    #
    # Inject outliers
    #
    indexOut <- matrix(nrow=j, ncol=m,
                       sample(c(-1,1,0), m*j, replace=TRUE, prob=c(freq/2, freq/2, 1-freq)))
    indexOut <- switch(match.arg(inj),
                       low  = -abs(indexOut),
                       high =  abs(indexOut),
                       indexOut
    )
    
    out_range <- c(0,1)
    n_rejected <- 0
    list_index <- which(indexOut != 0, arr.ind = TRUE)
    for(i in seq_len(nrow(list_index))){
      row <- list_index[i,'row']
      col <- list_index[i,'col']
      fc <- zScore * datasd[row,col]
      outlier_psi <- indexOut[row,col] * fc + lmu[row,col]
      art_out <- 2^outlier_psi
      
      if(art_out > out_range[1] && art_out < out_range[2]){
        # k/n = psi therefore k = psi * n (plus pseudocounts)
        if(type == "psi3"){
          k_new <- round( art_out * (n[row,col] + (2*pseudocount())) - pseudocount() )
          k[row,col] <- min(n[row,col], max(0, k_new) ) # max(0,...) to ensure k is never negative, min(n, ...) to ensure k <= n
        } else{
          k_new <- round( art_out * (nonSplit[row,col] + (2*pseudocount())) - pseudocount() )
          k[row,col] <- min(nonSplit[row,col], max(0, k_new) ) # max(0,...) to ensure k is never negative, min(n, ...) to ensure k <= n
        }
      }else{
        #remove outliers with psi < 0 or > 1
        indexOut[row,col] <- 0 
        n_rejected <- n_rejected + 1
      }
    }
    mode(k) <- "integer"
    
    setAssayMatrix(fds=fds, name="trueOutliers", type=type) <- indexOut
    counts(fds, type=type, side="ofInterest") <- k
    if(type == "psi3"){
      counts(fds, type=type, side="other") <- (n - k)
    } else{
      counts(fds, type=type, side="other") <- (nonSplit - k)
    }
    
    # for now: same values for psi3 and psi5
    if(type== "psi3"){
      counts(fds, type="psi5", side="ofInterest") <- k
      counts(fds, type="psi5", side="other") <- (n - k)
      setAssayMatrix(fds=fds, name="trueOutliers", type="psi5") <- indexOut
    }
  
  }
  
  return(fds)
}

# log2 fold change of measured psi vs AE predicted psi
log2fc <- function(realPsi, predictedPsi){
  return( log2(realPsi) - log2(predictedPsi) )
}
