
#' 
#' Create an simulated example data set for FraseR
#' 
#' Simulates a data set based on random counts following a 
#' negative binomial distribution with injected outliers with 
#' a fixed z score away from the mean of the gene.
#' 
#' @param j Number of simulated junctions 
#' @param m Number of simulated samples
#' @param q number of simulated latent variables.
#' 
#' @return An FraserDataSet containing an example dataset based on 
#'         simulated data
#' 
#' @examples
#' # A generic dataset 
#' fds1 <- makeSimulatedFraserDataSet()
#' fds1
#' fds1 <- injectOutliers(fds1, type="psi3", deltaPSI=0.4, freq=1E-3, inj="both")
#' 
#' # A generic dataset with specificed sample size and injection method
#' fds2 <- makeSimulatedFraserDataSet(m=100, j=500, inj='low')
#' fds2
#' 
#' @export
makeSimulatedFraserDataSet_BetaBinomial <- function(m=200, j=10000, q=10){
  
  # Simulate counts at each junction using negative binomial
  # nMean <- rnorm(j, 348, 10)                     # nbinom means for n (mean(N) on kremer dataset = 347.62, median=37)
  # theta <- rnorm(j, 0.27, 0.01)                  # nbinom dispersion (~0.27 on N from kremer dataset)
  
  # Simulate N = nonSplit + n
  junction_n <- rnbinom(j, mu=348, size=0.27)
  N <- t(vapply(junction_n, function(x){
                                n <- round(rnorm(n=m, sd=10, mean=x))
                                n <- pmax(n, 0)}, 
                double(m)))
  # N <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
  mode(N) <- 'integer'
  
  #
  # Simulate covariates for SE
  #
  sdVec <- rep(0.5, m)                           # sd for H matrix
  H_true <- matrix(rnorm(m*q, 0, 1), nrow=m, ncol=q)
  D_true <- matrix(rnorm(j*q, sd=sdVec), nrow=j, ncol=q)
  y_se <- D_true %*% t(cbind(H_true))
  mu_se     <- predictMuCpp(y_se)
  rho_se <- abs(rnorm(j, mean=0.0001, sd=0.05))     # betaBin dispersion
  
  # Draw nonSplit reads from BB
  nonSplit <-matrix(rbetabinom(j*m, size=N, prob=mu_se, rho=rho_se), nrow=j, ncol=m)
  mode(nonSplit) <- 'integer'
  
  # Set n = N - nonSplit
  n <- N - nonSplit
  
  #
  # Simulate covariates for PSI3=PSI5
  #
  sdVec <- rep(0.5, m)                           # sd for H matrix
  H_true <- matrix(rnorm(m*q, 0, 1), nrow=m, ncol=q)
  D_true <- matrix(rnorm(j*q, sd=sdVec), nrow=j, ncol=q)
  y_psi <- D_true %*% t(cbind(H_true))
  mu_psi     <- predictMuCpp(y_psi)
  rho_psi <- abs(rnorm(j, mean=0.0001, sd=0.05))     # betaBin dispersion
  
  # Draw nonSplit reads from BB
  k <-matrix(rbetabinom(j*m, size=n, prob=mu_psi, rho=rho_psi), nrow=j, ncol=m)
  mode(k) <- 'integer'
  
  
  #
  # Create FraseR data set
  #
  sampleIDs <- paste0("sample", seq_len(m))
  anno <- data.table(sampleID = sampleIDs, bamFile=rep(NA, m))
  fds <- FraseRDataSet(colData=anno)
  
  # put in n as rawcountsJ first so it doesn't complain later when assinging k to it
  junctionData <- SummarizedExperiment(
    colData=colData(fds),
    assays=list(rawCountsJ=k),
    rowRanges=GRanges(seqnames=rep("chr1", j), ranges=IRanges(start=seq_len(j), width=1 ))
  )
  nonSplitData <- SummarizedExperiment(
    colData=colData(fds),
    assays=list(rawCountsSS=nonSplit),
    rowRanges=GRanges(seqnames=rep("chr1", j), ranges=IRanges(start=seq_len(j), width=1 ))
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
  
  # Store "other" counts
  counts(fds, type="psi3", side="other") <-    n - k
  counts(fds, type="psi5", side="other") <-    n - k
  counts(fds, type="psiSite", side="other") <- n
  
  # store information about the simulation in the fds (same values for psi3 and psi5)
  for(type in c("psi3", "psi5")){
    setAssayMatrix(fds=fds, name="truePSI", type=type)      <- mu_psi 
    setAssayMatrix(fds=fds, name="trueLogitPSI", type=type) <- y_psi
    mcols(fds, type=type)[,"trueRho"]                       <- rho_psi
    
    # needed so that subsetting the fds works later
    mcols(fds, type=type)[["startID"]] <- 1:nrow(mcols(fds, type=type))
    mcols(fds, type=type)[["endID"]]   <- 1:nrow(mcols(fds, type=type))
  }
  
  # store info for SE 
  setAssayMatrix(fds=fds, name="truePSI", type="psiSite")      <- mu_se
  setAssayMatrix(fds=fds, name="trueLogitPSI", type="psiSite") <- y_se
  mcols(fds, type="psiSite")[,"trueRho"]                       <- rho_se
  
  # needed so that subsetting the fds works later
  mcols(fds, type="psiSite")[["startID"]]      <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["endID"]]        <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["spliceSiteID"]] <- 1:nrow(mcols(fds, type="psiSite"))
  
  return(fds)
  
}

#
# Generate a simulated fds using a Dirichlet-Mulitnomial distribution
#
makeSimulatedFraserDataSet_Multinomial <- function(m=200, j=1000, q=10){
  
  
  #
  # Simulate groups of junctions having the same donor/acceptor site
  #
  d <- round(j/4)              # nr of donors/acceptors
  donorGroups <- seq_len(d)    # ids of the groups (1,2,...,d)
  junctionGroups <- sort(sample(x = donorGroups, size=j, replace = TRUE))   # assign junction to donor/acceptor groups 
  
  # Find groups with only one junction and assign them randomly to another group
  singleDonors <- as.integer(names(which(table(junctionGroups) == 1)))
  usedGroups <- unique(junctionGroups)
  junctionGroups[which(junctionGroups %in% singleDonors)] <- sample(usedGroups[!(usedGroups %in% singleDonors)], size=length(singleDonors), replace = TRUE)
  junctionGroups <- sort(junctionGroups)
  
  # Set d to the actual number of groups (some might not have been assigned any junction at all)
  donorGroups <- donorGroups[donorGroups %in% unique(junctionGroups)]
  d <- length(donorGroups)
  
  # 
  # SE simulated as additional junction for every group
  #
  nrJunctions <- j                      # store for later
  psiJunctionGroups <- junctionGroups   # store for later
  j <- j + d                            # add one junction for every group
  junctionGroups <- sort(c(junctionGroups, donorGroups))
  
  
  
  #
  # Simulate n for each junction using negative binomial
  #
  # nMean <- rnorm(d, 348, 10)           # nbinom means for n (mean(N) on kremer dataset = 347.62, median=37)
  # theta <- rnorm(d, 0.27, 0.01)        # nbinom dispersion (~0.27 on N from kremer dataset)
  # n <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
  junction_n <- rnbinom(j, mu=348, size=0.27)
  n <- t(vapply(junction_n, function(x){
                              n <- round(rnorm(n=m, sd=10, mean=x))
                              n <- pmax(n, 0)}, 
                double(m)))
  
  # Sum up the values of n within one donor/acceptor group to get the n value of the group
  dt_n <- as.data.table(n)
  sum_n <- cbind(dt_n, junctionGroups)
  res <- sapply(colnames(dt_n),function(x){
    sum_n[,c(paste0(x, "_sum")):=sum(get(x)),by=junctionGroups]
  })
  
  # Extract final n matrix (all junctions within a group have the same n value)
  sum_cols <- paste0(colnames(dt_n), "_sum")
  n <- as.matrix(sum_n[,..sum_cols])
  colnames(n) <- NULL
  mode(n) <- 'integer'
  
  
  
  #
  # Simulate betaBin dispersion for every donor/acceptor group, same dispersion for all junctions within a group
  #
  group_rho <- abs(rnorm(d, mean=0.0001, sd=0.05))     # group dispersion
  dt1 <- data.table(group=donorGroups, rho=group_rho)
  dt2 <- data.table(group=junctionGroups)
  rhoFull <- merge(dt1, dt2, by="group")
  rho <- rhoFull[,rho]                           # set dispersion for every junction 
  
  
  
  #
  # Simulate covariates.
  #
  sdVec <- rep(0.5, m)                           # sd for H matrix
  H_true <- matrix(rnorm(m*q, 0, 1), nrow=m, ncol=q)
  D_true <- matrix(rnorm(j*q, sd=sdVec), nrow=j, ncol=q)
  y_true <- D_true %*% t(cbind(H_true))
  ae_mu     <- predictMuCpp(y_true)
  
  # Use softmax on mu to get mu within one group to sum up to 1
  softmax <- function(x){
    sum <- sum(exp(x))
    return( exp(x) / sum)
  }
  softmax_mu <- sapply(donorGroups, function(groupID){
    group_mu <- ae_mu[which(junctionGroups == groupID),]
    return( apply(group_mu, 2, softmax) )
  })
  mu <- do.call(rbind, softmax_mu)
  
  # Calcluate alpha values
  alpha <- mu * (1-rho)/rho
  
  
  
  # 
  # Simulate psi3=psi5 and SE jointly with Dirchilet-Multinomial
  # 
  res <- sapply(donorGroups, function(groupID){

    # Get the indices of the junctions in this group
    pos <- which(junctionGroups == groupID)
    
    # get alpha and n values for this group
    group_alpha <- alpha[pos,]
    group_n   <- n[pos,]
    
    # draw counts from a dirichlet multinomial distribution (jointly for PSi and SE)
    counts <- t(rdirmnom(m, group_n[1,], t(group_alpha)))
    
    # First row of resulting counts represents non split reads for SE, the other rows are the counts for the "real" junctions
    nonSplit <- counts[1,]
    k <- counts[-1,]
    
    # Substract nonSplit reads from total N to get n for PSI (=other for SE)
    group_n   <- group_n - matrix(nonSplit, nrow=nrow(group_n), ncol=ncol(group_n), byrow = TRUE)
    
    # Also get mu and rho for this group (to split it into PSI and SE part for storing later)
    group_mu  <- mu[pos,]
    group_rho <- rho[pos]
    
    # First row/value is the one for SE
    se_n     <- group_n[1,]
    se_mu    <- group_mu[1,]
    se_alpha <- group_alpha[1,]
    se_rho   <- group_rho[1]
    
    # The other rows/values are the ones relevant for PSI
    psi_n     <- group_n[-1,]
    psi_alpha <- group_alpha[-1,]
    psi_rho   <- group_rho[-1] 
    # scale mu's so that they sum up to 1 again (after mu for SE is removed)
    psi_mu    <- group_mu[-1,] /  matrix(colSums(group_mu[-1,]), nrow=nrow(group_mu[-1,]), ncol=ncol(group_mu[-1,]), byrow = TRUE)  
    
    # return everyting relevant for storing counts, n, alpha, mu, rho split into the parts for PSI and SE
    return(list(k=k, nonSplit=nonSplit, n=psi_n, group_n=se_n, psi_mu=psi_mu, se_mu=se_mu, 
                psi_alpha=psi_alpha, se_alpha=se_alpha, psi_rho=psi_rho, se_rho=se_rho))
  })
  
  # Extract k and nonSplit reads
  k <- do.call(rbind, res[1,])
  nonSplit <- do.call(rbind, res[2,])
  mode(k) <- 'integer'
  mode(nonSplit) <- 'integer'
  
  # Extract n value for each junction and group n value for SE
  psi_n <- do.call(rbind, res[3,])
  se_n <- do.call(rbind, res[4,])
  
  # Extract mu values for PSI and SE
  psi_mu <- do.call(rbind, res[5,])
  se_mu <- do.call(rbind, res[6,])
  
  # Extract alpha values for PSI and SE
  psi_alpha <- do.call(rbind, res[7,])
  se_alpha <- do.call(rbind, res[8,])
  
  # Extract rho values for PSI and SE
  psi_rho <- unlist(res[9,])
  se_rho <- unlist(res[10,])
  
  
  
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
    rowRanges=GRanges(seqnames=rep("chr1", nrJunctions), ranges=IRanges(start=psiJunctionGroups, width=1))
  )
  nonSplitData <- SummarizedExperiment(
    colData=colData(fds),
    assays=list(rawCountsSS=nonSplit),
    rowRanges=GRanges(seqnames=rep("chr1", d), ranges=IRanges(start=donorGroups, width=1))
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
  
  
  # store other counts (n-k) for psi3=psi5 and SE (other=n)
  counts(fds, type="psi3", side="other")    <- psi_n - k
  counts(fds, type="psi5", side="other")    <- psi_n - k
  counts(fds, type="psiSite", side="other") <- se_n
  
  for(type in c("psi3", "psi5")){
    # store information about the simulation in the fds (same values for psi3 and psi5)
    setAssayMatrix(fds=fds, name="truePSI", type=type)      <- psi_mu 
    setAssayMatrix(fds=fds, name="trueLogitPSI", type=type) <- qlogis(psi_mu) 
    setAssayMatrix(fds=fds, name="trueAlpha", type=type)    <- psi_alpha 
    mcols(fds, type=type)[,"trueRho"]                       <- psi_rho
    
    # needed so that subsetting the fds works later
    mcols(fds, type=type)[["startID"]] <- 1:nrow(mcols(fds, type=type))
    mcols(fds, type=type)[["endID"]]   <- 1:nrow(mcols(fds, type=type))
  }
  
  # store info for SE 
  setAssayMatrix(fds=fds, name="truePSI", type="psiSite")      <- se_mu
  setAssayMatrix(fds=fds, name="trueLogitPSI", type="psiSite") <- qlogis(se_mu) 
  setAssayMatrix(fds=fds, name="trueAlpha", type="psiSite")    <- se_alpha 
  mcols(fds, type="psiSite")[,"trueRho"]                       <- se_rho
  
  # needed so that subsetting the fds works later
  mcols(fds, type="psiSite")[["startID"]]      <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["endID"]]        <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["spliceSiteID"]] <- 1:nrow(mcols(fds, type="psiSite"))
  
  return(fds)
  
}

#
# Inject artificial outliers into a fds
#
injectOutliers <- function(fds, type=type, deltaPSI=0.3, freq=1E-3, inj=c('both', 'low', 'high'), method=c('samplePSI', 'meanPSI')){
  
  # Extract needed data from fds
  m <- ncol(fds)
  j <- nrow(mcols(fds, type=type))
  
  k <- K(fds, type=type)
  n <- N(fds, type=type)
  o <- counts(fds, type=type, side="other")
  
  psi <- (k + pseudocount())/(n + 2*pseudocount())
  
  if(type == "psi5"){
    o_psi5 <- o
    o_psi3 <- counts(fds, type="psi3", side="other")
  }
  else if(type =="psi3"){
    o_psi5 <- counts(fds, type="psi5", side="other")
    o_psi3 <- o
  }
  
  # Start and end positions of junctions
  dt <- data.table(id=1:nrow(fds),chr=as.character(seqnames(fds)), start=start(fds), end=end(fds))
  
  
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
  
  n_rejected <- 0
  list_index <- which(indexOut != 0, arr.ind = TRUE)
  for(i in seq_len(nrow(list_index))){
    row <- list_index[i,'row']
    col <- list_index[i,'col']
    n_ji <- n[row,col]
    o_ji <- o[row,col]
    
    #direction of outlier based on data psi value (low, high)
    indexOut[row,col] <- ifelse(psi[row,col] < 0.5, 1, -1)
    
    # new psi based on psi of sample i and junction j (=psi[j,i]) or based on junction mean = mean(psi[j,])
    outlier_psi <- switch(match.arg(method), 
                          samplePSI = psi[row,col] + indexOut[row,col] * deltaPSI,
                          meanPSI   = mean(psi[row,]) + indexOut[row,col] * deltaPSI)
    
    # change both k and n=k+o (and take pseudocounts into account): (k+1)/(k+o+2)=psi -> k = psi/(1-psi) * (o+2) - 1/(1-psi)
    k_new <- round( (outlier_psi*(o_ji + 2*pseudocount()) - pseudocount())/(1-outlier_psi) )
    # old way: change only k, n stays the same: k/n = psi therefore k = psi * n (plus pseudocounts)
    # k_new <- round( outlier_psi * (n_ji + (2*pseudocount())) - pseudocount() )
    
    if(outlier_psi > 1 || outlier_psi < 0 || # ensure that psi is between 0 and 1 
       k_new < 0 || k_new > n_ji){ # and ensure k is never negative and always <= n
      
      #remove outliers with psi < 0 or > 1 or k < 0 or k > n
      indexOut[row,col] <- 0 
      n_rejected <- n_rejected + 1
      
    }else{ # 0 <= psi <= 1 and 0 <= k <= n
      k[row,col] <- k_new  
      n_ji <- k_new + o_ji
      
      if(type == "psi5"){
        nPsi5 <- n_ji
        nPsi3 <- k_new + o_psi3[row,col]
      }
      else if(type == "psi3"){
        nPsi5 <- k_new + o_psi5[row,col]
        nPsi3 <- n_ji
      }
    
      # change other counts for all other junctions having the same start or end position
      if(type == "psi5" || type=="psi3"){
        sameDonor <- dt[start==dt[row,start] & chr==dt[row,chr], id]
        sameDonor <- sameDonor[sameDonor != row]
        
        # if there is at least one other junction starting at the same position
        if(length(sameDonor) > 0){
          for(junction in sameDonor){
            # adjust cother psi5 counts to correct value
            o_psi5[junction,col] <- nPsi5 - k[junction,col]
          }
        }
        
        sameAcceptor <- dt[end==dt[row,end] & chr==dt[row,chr], id]
        sameAcceptor <- sameAcceptor[sameAcceptor != row]
        
        # if there is at least one other junction ending at the same position
        if(length(sameAcceptor) > 0){
          for(junction in sameAcceptor){
            # adjust other psi3 counts to correct value
            o_psi3[junction,col] <- nPsi3 - k[junction,col]
          }
        }
      }
      
    }
    
  }
  
  # store positions of true outliers
  setAssayMatrix(fds=fds, name="trueOutliers", type=type) <- indexOut
  
  # store new k and o counts including the outlier counts (for SE: o doesn't change)
  counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- k
  if(type == "psi5"){
    counts(fds, type=type, side="other", HDF5=FALSE) <- o_psi5
    counts(fds, type="psi3", side="other", HDF5=FALSE) <- o_psi3
  }
  else if(type == "psi3"){   
    counts(fds, type=type, side="other", HDF5=FALSE) <- o_psi3
    counts(fds, type="psi5", side="other", HDF5=FALSE) <- o_psi5
  } 
  
  return(fds)
  
}


# log2 fold change of measured psi vs AE predicted psi
log2fc <- function(realPsi, predictedPsi){
  return( log2(realPsi) - log2(predictedPsi) )
}
