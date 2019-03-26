
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
makeSimulatedFraserDataSet_BetaBinomial <- function(m=200, j=10000, q=10, ...){

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
  fds <- FraseRDataSet(colData=anno, ...)

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
makeSimulatedFraserDataSet_Multinomial <- function(m=200, j=1000, q=10, ...){


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
  fds <- FraseRDataSet(colData=anno, ...)

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
    mcols(fds, type=type)[["startID"]] <- psiJunctionGroups # 1:nrow(mcols(fds, type=type))
    mcols(fds, type=type)[["endID"]]   <- psiJunctionGroups # 1:nrow(mcols(fds, type=type))
  }

  # store info for SE
  setAssayMatrix(fds=fds, name="truePSI", type="psiSite")      <- se_mu
  setAssayMatrix(fds=fds, name="trueLogitPSI", type="psiSite") <- qlogis(se_mu)
  setAssayMatrix(fds=fds, name="trueAlpha", type="psiSite")    <- se_alpha
  mcols(fds, type="psiSite")[,"trueRho"]                       <- se_rho

  # needed so that subsetting the fds works later
  mcols(fds, type="psiSite")[["startID"]]      <- donorGroups # 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["endID"]]        <- 1:nrow(mcols(fds, type="psiSite"))
  mcols(fds, type="psiSite")[["spliceSiteID"]] <- donorGroups # 1:nrow(mcols(fds, type="psiSite"))

  return(fds)

}


#
# Inject artificial outliers in an existing fds
#
injectOutliers <- function(fds, type=type, nrOutliers=500, deltaPSI=pmin(0.7, pmax(-0.7, rnorm(nrOutliers, 0, 0.5))),
                           method=c('meanPSI', 'samplePSI', 'simulatedPSI'), swap=TRUE, verbose=FALSE){

  # copy original k and o
  if(type == "psiSite"){
    setAssayMatrix(fds, type="psiSite", "originalCounts") <- counts(fds, type="psiSite", side="ofInterest")
    setAssayMatrix(fds, type="psiSite", "originalOtherCounts") <- counts(fds, type="psiSite", side="other")
  }
  else{
    setAssayMatrix(fds, type="psi5", "originalCounts") <- counts(fds, type="psi5", side="ofInterest")
    setAssayMatrix(fds, type="psi5", "originalOtherCounts") <- counts(fds, type="psi5", side="other")
    setAssayMatrix(fds, type="psi3", "originalOtherCounts") <- counts(fds, type="psi3", side="other")
  }

  # get infos from the fds
  m <- ncol(fds)
  j <- nrow(mcols(fds, type=type))

  k <- as.matrix(K(fds, type=type))
  n <- as.matrix(N(fds, type=type)) # as.matrix(..) needed so that n doesn't change after new k is stored (needed only for swap=FALSE)
  o <- as.matrix(counts(fds, type=type, side="other"))

  psi <- (k + pseudocount())/(n + 2*pseudocount())

  # nedded to get the junction counts when type is SE
  if(type == "psiSite"){
    k_other <- K(fds, type = "psi5")
  }

  # get psi from the simulation if simulated data is used
  if(match.arg(method) == 'simulatedPSI'){
    truePSI <- getAssayMatrix(fds, "truePSI", type=type)
  }

  # matrices to store indices of outliers and their delta psi
  indexOut <- matrix(0, nrow = j, ncol = m)
  indexDeltaPSI <- matrix(0, nrow = j, ncol = m)

  # junctions to choose randomly from
  available_junctions <- seq_len(j)

  # loop over all outliers that should be injected
  for(i in seq_len(nrOutliers)){

    successful <- FALSE

    # for each outlier, draw random junction-sample pair and check if outlier can be injected there
    while(!successful){

      junction <- sample(available_junctions, 1)
      sample <- sample(m, 1)

      n_ji <- n[junction,sample]

      if(n_ji < 10){ next }

      # new psi based on psi of sample i and junction j (=psi[j,i]) or based on junction mean = mean(psi[j,]) or based on the psi used during the simulation
      outlier_psi <- switch(match.arg(method),
                            samplePSI    = psi[junction,sample] + deltaPSI[i],
                            meanPSI      = mean(psi[junction,]) + deltaPSI[i],
                            simulatedPSI = truePSI[junction,sample] + deltaPSI[i])

      # change k based on n and outlier_psi (and take pseudocounts into account): (k+1)/(n+2)=psi -> k = psi*(n+2) - 1
      k_new <- round( outlier_psi*(n_ji + 2*pseudocount()) - pseudocount() )

      # try again for outliers with psi < 0 or > 1 or k < 0 or k > n
      if(outlier_psi <= 1 && outlier_psi >= 0 && k_new >= 0 && k_new <= n_ji){

        if(swap){
          #
          # check if junction can be swapped with other junction from same donor/acceptor
          #
          if(type == "psi5"){
            samePos <- which(mcols(fds, type=type)[["startID"]] == mcols(fds, type=type)[["startID"]][junction])
          }
          else if(type == "psi3"){
            samePos <- which(mcols(fds, type=type)[["endID"]] == mcols(fds, type=type)[["endID"]][junction])
          }
          else{
            # for SE: swap with one junction at the same position
            samePos <- which(mcols(fds, type="psi5")[["startID"]] == mcols(fds, type="psiSite")[["startID"]][junction])
            samePos <- c(samePos, which(mcols(fds, type="psi3")[["endID"]] == mcols(fds, type="psiSite")[["endID"]][junction]))
            # n for other junctions is o for SE
            n_ji <- o[junction, sample]
          }
          group_others <- samePos[samePos != junction]

          if(length(group_others) == 0){
            # no other junction to swap with
            next
          }

          k_group <- switch(type,
                            "psiSite" = k_other[group_others, sample],
                            k[group_others, sample])

          # difference between "old" and "new" k
          diff <- k[junction, sample] - k_new

          # try if swapping counts with other junction at this position is possible
          g <- 1
          foundSwap <- FALSE
          while(!foundSwap){

            k_group_new <- k_group[g] + diff

            if(k_group_new >= 0 && k_group_new <= n_ji){

              if(type == "psiSite"){
                k_other[group_others[g],sample] <- k_group_new
              }
              else{
                k[group_others[g],sample] <- k_group_new
              }

              k[junction, sample] <- k_new

              foundSwap <- TRUE
              successful <- TRUE

              # also store position of other junction used in swap
              indexOut[group_others[g], sample] <- ifelse(deltaPSI[i] >= 0, -2, 2)
              indexDeltaPSI[junction, sample] <- deltaPSI[i]

            }

            g <- g + 1
            if(g > length(group_others)){
              break
            }
          }

        }
        else{ # no swap should be done (e.g for BB simulated data)

          # just set k to computed k_new
          k[junction, sample] <- k_new
          successful <- TRUE
        }

      }

      if(successful){
        indexOut[junction, sample] <- ifelse(deltaPSI[i] >= 0, 1, -1)
        indexDeltaPSI[junction, sample] <- deltaPSI[i]
      }

    }

    # remove junction where outlier was injected from available junctions (each junction at most one outlier)
    available_junctions <- available_junctions[available_junctions != junction]

    if(verbose){
      print(paste("Injected outlier", i, "at index", junction, ",", sample, ": new k =", k_new, ", n =", n_ji))
    }

  }

  # store positions of true outliers and their true delta PSI value
  setAssayMatrix(fds=fds, name="trueOutliers", type=type) <- indexOut
  metadata(fds)[[paste0("trueDeltaPSI_",type)]] <- indexDeltaPSI

  # store new k counts which include the outlier counts
  counts(fds, type=type, side="ofInterest") <- k

  # for SE: also store modified junction counts
  if(type == "psiSite"){
    counts(fds, type="psi5", side="ofInterest") <- k_other
  }

  # re-calculate other counts if swapping was done
  if(swap){
    fds <- calculatePSIValues(fds, overwriteCts = TRUE)
  }
  else{
    # when no swapping was done, just put other counts as n-k
    counts(fds, type=type, side="other") <- n-k
  }

  return(fds)

}


# log2 fold change of measured psi vs AE predicted psi
log2fc <- function(realPsi, predictedPsi){
  return( log2(realPsi) - log2(predictedPsi) )
}
