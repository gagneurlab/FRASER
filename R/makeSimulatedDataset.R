
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
#' @param deltaPSI delta PSI change of in-silico outliers.
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
makeSimulatedFraserDataSet <- function(m=200, j=10000, q=10, freq=1E-3, deltaPSI=0.2,
                                       inj=c('both', 'low', 'high')){
  
  # nMean <- ceiling(rlnorm(j, 5)) + 10     
  # theta <- rlnorm(j, 4)                         
  # nMean <- rnbinom(j, 348, 0.27)               # nbinom means for n (mean(N) on kremer dataset = 347.62, median=37)
  nMean <- rnorm(j, 348, 10)                     # nbinom means for n (mean(N) on kremer dataset = 347.62, median=37)
  theta <- rnorm(j, 0.27, 0.01)                  # nbinom dispersion (~0.27 on N from kremer dataset)
  
  # simulate counts at each junction
  n <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
  # nonSplit <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
  # n <- matrix(abs(round(rnorm(j*m, mean=nMean, sd=nMean/4))), nrow=j, ncol=m)
  # nonSplit <- matrix(abs(round(rnorm(j*m, mean=nMean, sd=nMean/4))), nrow=j, ncol=m)
  
  
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
    assays=list(rawCountsSS=n),
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
    H_true <- matrix(rnorm(m*q, 0, 2), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, sd=sdVec), nrow=j, ncol=q)
    y_true <- D_true %*% t(cbind(H_true))
    mu     <- predictMuCpp(y_true)
  
    # Simulate counts
    if(type== "psi3"){
      k <-matrix(rbetabinom(j*m, size=n, prob=mu, rho=rho), nrow=j, ncol=m)
    } else{
      k <-matrix(rbetabinom(j*m, size=n, prob=mu, rho=rho), nrow=j, ncol=m)
      # mu <- mu/2
      # y_true <- qlogis(mu)
      # k <-matrix(rbetabinom(j*m, size=n, prob=(mu/(1-mu)), rho=rho), nrow=j, ncol=m)
    }
    mode(n) <- 'integer'
    mode(k) <- 'integer'
  
    psi <- (k + pseudocount())/(n + 2*pseudocount())
    # l2fc <- log2fc(psi, mu)
    # datasd <- rowSds(l2fc)
    # datasd <- matrix(datasd, nrow=j, ncol=m)
    # lmu <- log2(mu)
    
    setAssayMatrix(fds=fds, name="truePSI", type=type) <- mu 
    setAssayMatrix(fds=fds, name="trueLogitPSI", type=type) <- y_true 
    # setAssayMatrix(fds=fds, name="trueSd", type=type)  <- datasd
    mcols(fds, type=type)[,"trueRho"] <- rho
    
    # needed so that subsetting the fds works later
    mcols(fds, type=type)[["startID"]] <- 1:nrow(mcols(fds, type=type))
    mcols(fds, type=type)[["endID"]] <- 1:nrow(mcols(fds, type=type))
    
    # for now: same values for psi3 and psi5
    if(type== "psi3"){
      setAssayMatrix(fds=fds, name="truePSI", type="psi5") <- mu 
      setAssayMatrix(fds=fds, name="trueLogitPSI", type="psi5") <- y_true
      # setAssayMatrix(fds=fds, name="trueSd", type="psi5")  <- datasd
      mcols(fds, type="psi5")[,"trueRho"] <- rho
      
      # needed so that subsetting the fds works later
      mcols(fds, type="psi5")[["startID"]] <- 1:nrow(mcols(fds, type=type))
      mcols(fds, type="psi5")[["endID"]] <- 1:nrow(mcols(fds, type=type))
      
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
      # fc <- zScore * datasd[row,col]
      # outlier_psi <- indexOut[row,col] * fc + lmu[row,col]
      # art_out <- 2^outlier_psi
      art_out <- mu[row,col] + indexOut[row,col] * deltaPSI  # mu is true simulated psi
      
      # k/n = psi therefore k = psi * n (plus pseudocounts)
      if(type == "psi3"){
        n_pos <- n[row,col]
      } else{
        n_pos <- nonSplit[row,col]
      }
      k_new <- round( art_out * (n_pos + (2*pseudocount())) - pseudocount() )
      
      if(art_out > out_range[1] && art_out < out_range[2] && 
         k_new >= 0 && k_new <= n_pos){ # max(0,...) to ensure k is never negative, min(n, ...) to ensure k <= n
          
        k[row,col] <- k_new  
        
      }else{
        #remove outliers with psi < 0 or > 1 or k < 0 or k > n
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
      counts(fds, type=type, side="other") <- n - k # just n? (nonSplit=k)
    }
    
    # for now: same values for psi3 and psi5
    if(type== "psi3"){
      # counts(fds, type="psi5", side="ofInterest") <- k
      counts(fds, type="psi5", side="other") <- (n - k)
      setAssayMatrix(fds=fds, name="trueOutliers", type="psi5") <- indexOut
    }
  
  }
  
  return(fds)
}

# Generate a simulated fds using a Dirichlet-Mulitnomial distribution
makeSimulatedFraserDataSet_Multinomial <- function(m=200, j=10000, q=10){
  
  # nr of donors/acceptor
  d <- round(j/2)  
  donorGroups <- seq_len(d)
  # assign junction to donor/acceptor groups (making sure each donor/accptor has at least one junction)
  junctionGroups <- sort(c(donorGroups, sample(x = donorGroups, size=j-d, replace = TRUE)))
  
  # simulate n for each donor/acceptor group
  nMean <- rnorm(d, 348, 10)                     # nbinom means for n (mean(N) on kremer dataset = 347.62, median=37)
  theta <- rnorm(d, 0.27, 0.01)                  # nbinom dispersion (~0.27 on N from kremer dataset)
  group_n <- matrix(rnbinom(d*m, mu=nMean, size=theta), nrow=d, ncol=m)
  
  # construct full n matrix (all junctions within a donor/acceptor group have the same n)
  n <- group_n[rep(donorGroups, table(junctionGroups)),]
  mode(n) <- 'integer'
  
  #
  # Simulate covariates.
  #
  sdVec <- rep(0.5, m)                           # sd for H matrix
  H_true <- matrix(rnorm(m*q, 0, 1), nrow=m, ncol=q)
  D_true <- matrix(rnorm(j*q, sd=sdVec), nrow=j, ncol=q)
  y_true <- D_true %*% t(cbind(H_true))
  mu     <- predictMuCpp(y_true)
  
  rho <- abs(rnorm(j, mean=0.0001, sd=0.05))     # betaBin dispersion
  alpha <- mu * (1-rho)/rho
  

  # simulate psi 3 (=psi5 for now)
  res <- sapply(seq_len(d), function(groupID){
    
            group_alpha <- alpha[which(junctionGroups == groupID),]
            if(!is.matrix(group_alpha)){
              group_alpha <- matrix(group_alpha, nrow=1, ncol=length(group_alpha))
            }
            # draw vector of probablities for the multinomial distribution from a dirchlet distribution
            p <- VGAM::rdiric(1, t(group_alpha), is.matrix=TRUE)
            p[p == 0] <- 0.5 # p=0 causes problems in rmultinom, need to fix that somehow
            
            # draw k's from a multinomial distribution using each sample's n value
            k <- vapply(seq_len(m), function(i){ rmultinom(1, group_n[groupID,i], p[i,]) }, integer(ncol(p)))
            
            p[p == 1] <- 1 - runif(sum(p == 1), 1e-20, 1e-2) # p=1 causes problems with logit transformation, so p is set to some random value close to 1 for now
            
            return(list(k=k,p=t(p)))
  })
 
  new_mu <- do.call(rbind, res[2,])    # probablities used when drawing from the multinomial
  k <- do.call(rbind, res[1,])         # counts generated by drawing from the mulitnomial
  mode(k) <- 'integer'
  
  # # calculate psi3 of simulated counts
  # psi <- (k + pseudocount())/(n + 2*pseudocount())
  
  # simulate SE count (= nonSplit reads) (how to do this right?)
  mu_se <- mu
  rho_se <- abs(rnorm(j, mean=0.0001, sd=0.05))     # betaBin dispersion
  nonSplit <-matrix(rbetabinom(j*m, size=n, prob=mu_se, rho=rho_se), nrow=j, ncol=m)
  mode(nonSplit) <- 'integer'
  
  #
  # Create FraseR data set
  #
  sampleIDs <- paste0("sample", seq_len(m))
  anno <- data.table(sampleID = sampleIDs, bamFile=rep(NA, m))
  fds <- FraseRDataSet(colData=anno)
  
  # put in k as rawcountsJ (= nr reads spanning each junction)
  junctionData <- SummarizedExperiment(
    colData=colData(fds),
    assays=list(rawCountsJ=k),
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
  
  
  # store other counts (n-k) for psi3 (=psi5 for now) and SE (n - nonSplit)
  counts(fds, type="psi3", side="other") <- (n - k)
  counts(fds, type="psi5", side="other") <- (n - k)
  counts(fds, type="psiSite", side="other") <- n - nonSplit # just n?
  
  # store information about the simulation in the fds
  setAssayMatrix(fds=fds, name="truePSI", type="psi3") <- new_mu 
  setAssayMatrix(fds=fds, name="trueLogitPSI", type="psi3") <- qlogis(new_mu) 
  mcols(fds, type="psi3")[,"trueRho"] <- rho
  
  # needed so that subsetting the fds works later
  mcols(fds, type="psi3")[["startID"]] <- 1:nrow(mcols(fds, type="psi3"))
  mcols(fds, type="psi3")[["endID"]] <- 1:nrow(mcols(fds, type="psi3"))
  
  # for now: same values for psi3 and psi5
  setAssayMatrix(fds=fds, name="truePSI", type="psi5") <- new_mu 
  setAssayMatrix(fds=fds, name="trueLogitPSI", type="psi5") <- qlogis(new_mu) 
  mcols(fds, type="psi5")[,"trueRho"] <- rho
  mcols(fds, type="psi5")[["startID"]] <- 1:nrow(mcols(fds, type="psi5"))
  mcols(fds, type="psi5")[["endID"]] <- 1:nrow(mcols(fds, type="psi5"))
  
  # store info for SE 
  setAssayMatrix(fds=fds, name="truePSI", type="psiSite") <- mu_se 
  setAssayMatrix(fds=fds, name="trueLogitPSI", type="psiSite") <- qlogis(mu_se) 
  mcols(fds, type="psiSite")[,"trueRho"] <- rho_se
  
  # needed so that subsetting the fds works later
  mcols(fds, type="psiSite")[["startID"]] <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["endID"]] <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["spliceSiteID"]] <- 1:nrow(mcols(fds, type="psiSite"))
  
  return(fds)
  
}

# log2 fold change of measured psi vs AE predicted psi
log2fc <- function(realPsi, predictedPsi){
  return( log2(realPsi) - log2(predictedPsi) )
}
