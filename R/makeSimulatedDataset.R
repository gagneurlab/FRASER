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
#' @param distribution Either "BB" for a beta-binomial simulation or "DM" for a 
#' dirichlet-multinomial simulation.
#' @param ... Further arguments used to construct the FraseRDataSet.
#'
#' @return An FraserDataSet containing an example dataset based on
#'         simulated data
#'
#' @examples
#' # A generic dataset
#' fds1 <- makeSimulatedFraserDataSet()
#' fds1
#' fds1 <- injectOutliers(fds1, minDpsi=0.2, freq=1E-3)
#'
#' # A generic dataset with specificed sample size and injection method
#' fds2 <- makeSimulatedFraserDataSet(m=100, j=500, q=5)
#' fds2
#'
#' @rdname makeSimulatedFraserDataSet
#' @aliases makeSimulatedFraserDataSet, injectOutliers
#' @export
makeSimulatedFraserDataSet <- function(m=100, j=500, q=10,
                    distribution=c("BB", "DM"), ...){
    distribution <- match.arg(distribution)
    fds <- switch (distribution,
        BB = makeSimulatedFraserDataSet_BetaBinomial(m=m, j=j, q=q, ...),
        DM = makeSimulatedFraserDataSet_Multinomial(m=m, j=j, q=q, ...))

    fds
}

makeSimulatedFraserDataSet_BetaBinomial <- function(m=200, j=10000, q=10,
                    nbMeanLog=8, nbSize=0.55, rhoMeanlog=-5, rhoSdlog=1,
                    Hmean=0, Hsd=100, Dmean=0, Dsd=0.007, ...){

    # Simulate total coverage N = nonSplit + n at each junction using lognormal
    junction_n <- rlnorm(j, meanlog=nbMeanLog, sdlog=nbSize)
    N <- t(vapply(junction_n, function(x){
        # draw counts for every sample from negative binomial
        n <- rnbinom(m, mu=x, size=nbSize)},
        double(m)))
    # N <- matrix(rnbinom(j*m, mu=nbMean, size=nbSize), nrow=j, ncol=m)
    mode(N) <- 'integer'

    #
    # Simulate covariates for SE
    #
    H_true <- matrix(rnorm(m*q, mean=Hmean, sd=Hsd), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, mean=Dmean, sd=Dsd), nrow=j, ncol=q)
    b_true <- sample(c(rnorm(round(j*(1/6)), mean=3.5, sd=1.5),
            rnorm(round(j*(5/6)), mean=-2.5, sd=2.5)))
    y_se   <- t(predictYCpp(H_true, D_true, b_true))
    mu_se     <- predictMuCpp(y_se)
    # betaBin dispersion
    rho_se <- rlnorm(j, meanlog=rhoMeanlog, sdlog=rhoSdlog)

    # Draw nonSplit reads from BB
    nonSplit <-matrix(rbetabinom(j*m, size=N, prob=mu_se, rho=rho_se),
            nrow=j, ncol=m)
    mode(nonSplit) <- 'integer'

    # Set n = N - nonSplit
    n <- N - nonSplit

    #
    # Simulate covariates for PSI3=PSI5
    #
    H_true <- matrix(rnorm(m*q, mean=Hmean, sd=Hsd), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, mean=Dmean, sd=Dsd), nrow=j, ncol=q)
    b_true <- sample(c(rnorm(round(j*(1/6)), mean=-2.5, sd=2.5),
            rnorm(round(j*(5/6)), mean=3.5, sd=1.5)))
    y_psi  <- t(predictYCpp(H_true, D_true, b_true))
    mu_psi     <- predictMuCpp(y_psi)
    # betaBin dispersion
    rho_psi <- rlnorm(j, meanlog=rhoMeanlog, sdlog=rhoSdlog)

    # Draw nonSplit reads from BB
    k <-matrix(rbetabinom(j*m, size=n, prob=mu_psi, rho=rho_psi),
            nrow=j, ncol=m)
    mode(k) <- 'integer'


    #
    # Create FraseR data set
    #
    sampleIDs <- paste0("sample", seq_len(m))
    anno <- data.table(sampleID = sampleIDs, bamFile=rep(NA, m))
    fds <- FraseRDataSet(colData=anno, ...)

    # put in n as rawcountsJ first so it doesn't complain later
    # when assinging k to it
    junctionData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsJ=k),
        rowRanges=GRanges(seqnames=rep("chr1", j),
                ranges=IRanges(start=seq_len(j), width=1 ))
    )
    nonSplitData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsSS=nonSplit),
        rowRanges=GRanges(seqnames=rep("chr1", j),
                ranges=IRanges(start=seq_len(j), width=1 ))
    )
    fds <- new("FraseRDataSet",
            junctionData,
            name            = name(fds),
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

    # store information about the simulation in the fds
    # (same values for psi3 and psi5)
    for(type in c("psi3", "psi5")){
        setAssayMatrix(fds=fds, name="truePSI", type=type)      <- mu_psi
        setAssayMatrix(fds=fds, name="trueLogitPSI", type=type) <- y_psi
        mcols(fds, type=type)[,"trueRho"]                       <- rho_psi

        # needed so that subsetting the fds works later
        mcols(fds, type=type)[["startID"]] <- 
            seq_len(nrow(mcols(fds, type=type)))
        mcols(fds, type=type)[["endID"]]   <- 
            seq_len(nrow(mcols(fds, type=type)))
    }

    # store info for SE
    setAssayMatrix(fds=fds, name="truePSI", type="psiSite")      <- mu_se
    setAssayMatrix(fds=fds, name="trueLogitPSI", type="psiSite") <- y_se
    mcols(fds, type="psiSite")[,"trueRho"]                       <- rho_se

    # needed so that subsetting the fds works later
    siteIDs <- seq_row(mcols(fds, type="psiSite"))
    mcols(fds, type="psiSite")[["startID"]]      <- siteIDs
    mcols(fds, type="psiSite")[["endID"]]        <- siteIDs
    mcols(fds, type="psiSite")[["spliceSiteID"]] <- siteIDs

    return(fds)

}

#
# Generate a simulated fds using a Dirichlet-Mulitnomial distribution
#
makeSimulatedFraserDataSet_Multinomial <- function(m=200, j=1000, q=10,
                    groups=round(j*0.65), nbMeanLog=8, nbSize=0.55,
                    rhoMeanlog=-5, rhoSdlog=1, Hmean=0, Hsd=100, Dmean=0,
                    Dsd=0.007, ...){

    #
    # Simulate groups of junctions having the same donor/acceptor site
    #
    d <- groups                  # nr of donors/acceptors
    donorGroups <- seq_len(d)    # ids of the groups (1,2,...,d)
    # assign junction to donor/acceptor groups
    junctionGroups <- sort(sample(x = donorGroups, size=j, replace = TRUE))   

    # Set d to the actual number of groups (some might not have been assigned 
    # any junction at all)
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
    # Simulate n for each junction using lognormal and negative binomial
    #
    junction_n <- rlnorm(j, meanlog=nbMeanLog, sdlog=nbSize)
    n <- t(vapply(junction_n, function(x){
        # draw counts for every sample from negative binomial
        n <- rnbinom(m, mu=x, size=nbSize)},
        double(m)))
    # n <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)

    # Sum up the values of n within one donor/acceptor group to get the n value 
    # of the group
    dt_n <- as.data.table(n)
    sum_n <- cbind(dt_n, junctionGroups)
    res <- lapply(colnames(dt_n),function(x){
        sum_n[,c(paste0(x, "_sum")):=sum(get(x)),by=junctionGroups]
    })

    # Extract final n matrix (all junctions within a group have the same n 
    # value)
    sum_cols <- paste0(colnames(dt_n), "_sum")
    n <- as.matrix(sum_n[,..sum_cols])
    colnames(n) <- NULL
    mode(n) <- 'integer'



    #
    # Simulate betaBin dispersion for every donor/acceptor group, same 
    # dispersion for all junctions within a group
    #
    group_rho <- rlnorm(d, meanlog=rhoMeanlog, sdlog=rhoSdlog)# group dispersion
    # group_rho <- abs(rnorm(d, mean=0.0001, sd=0.05))     # group dispersion
    dt1 <- data.table(group=donorGroups, rho=group_rho)
    dt2 <- data.table(group=junctionGroups)
    rhoFull <- merge(dt1, dt2, by="group")
    rho <- rhoFull[,rho]                    # set dispersion for every junction



    #
    # Simulate covariates.
    #
    H_true <- matrix(rnorm(m*q, mean=Hmean, sd=Hsd), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, mean=Dmean, sd=Dsd), nrow=j, ncol=q)
    b_true <- double(j)
    y_true <- t(predictYCpp(H_true, D_true, b_true))
    ae_mu     <- predictMuCpp(y_true)

    # Use softmax on mu to get mu within one group to sum up to 1
    softmax <- function(x){
        sum <- sum(exp(x))
        return( exp(x) / sum)
    }
    softmax_mu <- lapply(donorGroups, function(groupID){
        group_mu <- ae_mu[which(junctionGroups == groupID),]
        return( apply(group_mu, 2, softmax) )
    })
    mu <- do.call(rbind, softmax_mu)

    # Calcluate alpha values
    alpha <- mu * (1-rho)/rho



    #
    # Simulate psi3=psi5 and SE jointly with Dirchilet-Multinomial
    #
    res <- vapply(donorGroups, function(groupID){

        # Get the indices of the junctions in this group
        pos <- which(junctionGroups == groupID)

        # get alpha and n values for this group
        group_alpha <- alpha[pos,]
        group_n   <- n[pos,]

        # draw counts from a dirichlet multinomial distribution 
        # (jointly for PSi and SE)
        counts <- t(rdirmnom(m, group_n[1,], t(group_alpha)))

        # First row of resulting counts represents non split reads for SE, 
        # the other rows are the counts for the "real" junctions
        nonSplit <- counts[1,]
        k <- counts[-1,]

        # Substract nonSplit reads from total N to get n for PSI (=other for SE)
        group_n   <- group_n - matrix(nonSplit, nrow=nrow(group_n), 
                                        ncol=ncol(group_n), byrow = TRUE)

        # Also get mu and rho for this group (to split it into PSI and SE part 
        # for storing later)
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
        psi_mu <- group_mu[-1,]
        if(is.null(nrow(psi_mu))){  # only one junction in the group, convert 
                                    # to matrix with one row so that colSums 
                                    # works
            psi_mu <- matrix(psi_mu, nrow=1)
        }
        # scale mu's so that they sum up to 1 again (after mu for SE is removed)
        psi_mu    <- psi_mu /  matrix(colSums(psi_mu), nrow=nrow(psi_mu), 
                                        ncol=ncol(psi_mu), byrow = TRUE)

        # return everyting relevant for storing counts, n, alpha, mu, 
        # rho split into the parts for PSI and SE
        return(list(k=k, nonSplit=nonSplit, n=psi_n, group_n=se_n, 
                    psi_mu=psi_mu, se_mu=se_mu, psi_alpha=psi_alpha, 
                    se_alpha=se_alpha, psi_rho=psi_rho, se_rho=se_rho))
    }, FUN.VALUE=list(k=numeric(m), nonSplit=numeric(m), n=numeric(m), 
                        group_n=numeric(m), psi_mu=numeric(m), se_mu=numeric(m), 
                        psi_alpha=numeric(m), se_alpha=numeric(m), 
                        psi_rho=numeric(1), se_rho=numeric(1)))

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
        rowRanges=GRanges(seqnames=rep("chr1", nrJunctions), 
                            ranges=IRanges(start=psiJunctionGroups, width=1))
    )
    nonSplitData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsSS=nonSplit),
        rowRanges=GRanges(seqnames=rep("chr1", d), 
                            ranges=IRanges(start=donorGroups, width=1))
    )
    fds <- new("FraseRDataSet",
                junctionData,
                name            = name(fds),
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
        # store information about the simulation in the fds 
        # (same values for psi3 and psi5)
        setAssayMatrix(fds=fds, name="truePSI", type=type)      <- psi_mu
        setAssayMatrix(fds=fds, name="trueLogitPSI", type=type) <- 
            qlogis(psi_mu)
        setAssayMatrix(fds=fds, name="trueAlpha", type=type)    <- psi_alpha
        mcols(fds, type=type)[,"trueRho"]                       <- psi_rho

        # needed so that subsetting the fds works later
        mcols(fds, type=type)[["startID"]] <- 
            psiJunctionGroups # 1:nrow(mcols(fds, type=type))
        mcols(fds, type=type)[["endID"]]   <- 
            psiJunctionGroups # 1:nrow(mcols(fds, type=type))
    }

    # store info for SE
    setAssayMatrix(fds=fds, name="truePSI", type="psiSite")      <- se_mu
    setAssayMatrix(fds=fds, name="trueLogitPSI", type="psiSite") <- 
        qlogis(se_mu)
    setAssayMatrix(fds=fds, name="trueAlpha", type="psiSite")    <- se_alpha
    mcols(fds, type="psiSite")[,"trueRho"]                       <- se_rho

    # needed so that subsetting the fds works later
    mcols(fds, type="psiSite")[["startID"]]      <- 
        donorGroups # 1:nrow(mcols(fds, type="psiSite"))
    mcols(fds, type="psiSite")[["endID"]]        <- 
        seq_len(nrow(mcols(fds, type="psiSite")))
    mcols(fds, type="psiSite")[["spliceSiteID"]] <- 
        donorGroups # 1:nrow(mcols(fds, type="psiSite"))

    return(fds)

}

#'
#' Inject artificial outliers in an existing fds
#' 
#' @param fds FraseRDataSet
#' @param type The psi type
#' @param freq The injection frequency.
#' @param minDpsi The minimal delta psi with which outliers will be injected.
#' @param minCoverage The minimal total coverage (i.e. N) required for a 
#' junction to be considered for injection of an outlier.
#' @param deltaDistribution The distribution from which the delta psi value of 
#' the injections is drawn (default: uniform distribution).
#' @param verbose Should additional information be printed during computation?
#' @param method Defines by which method the new psi of injections is computed, 
#' i.e. to which value the delta psi of the injection is added: "meanPSI" for 
#' adding to the mean psi of the junction over all samples or "samplePSI" to add 
#' to the psi value of the junction in the specific sample. "simulatedPSI" is 
#' only possible if a simulated dataset is used.
#' 
#' @return FraseRDataSet
#' @export
injectOutliers <- function(fds, type=c("psi5", "psi3", "psiSite"),
                    freq=1E-3, minDpsi=0.2, minCoverage=2,
                    deltaDistr="uniformDistr", verbose=FALSE,
                    method=c('meanPSI', 'samplePSI', 'simulatedPSI')){
    type <- match.arg(type, several.ok=TRUE)
    method <- match.arg(method)
    
    if(length(type) > 1){
        for(t in type){
            fds <- injectOutliers(fds, type=t, freq=freq, minDpsi=minDpsi,
                    minCoverage=minCoverage, deltaDistr=deltaDistr,
                    verbose=verbose, method=method)
        }
        return(fds)
    }
    
    # get infos from the fds
    if(paste0("originalCounts_", type) %in% assayNames(fds)){
        message("Use existing original counts and reinject outliers.")
        k <- as.matrix(getAssayMatrix(fds, type=type, "originalCounts"))
        o <- as.matrix(getAssayMatrix(fds, type=type, "originalOtherCounts"))
    } else {
        k <- as.matrix(counts(fds, type=type, side="ofI"))
        o <- as.matrix(counts(fds, type=type, side="other"))
    }
    n <- k + o
    m <- ncol(k)
    j <- nrow(k)
    
    # compute psi
    psi <- switch(method,
            samplePSI    = (k + pseudocount())/(n + 2*pseudocount()),
            meanPSI      = matrix(nrow=j, ncol=m,
                    rowMeans( (k + pseudocount())/(n + 2*pseudocount()) )),
            simulatedPSI = getAssayMatrix(fds, "truePSI", type=type))
    
    # 
    # matrix of donor/acceptor for possible injection
    # 
    tmpIndex <- getSiteIndex(fds, type)
    dt <- as.data.table(rowRanges(fds, type=type))[,.(
            chr=seqnames, start, end, strand, idxInCount=.I, 
            junctionID=tmpIndex)]
    dt[,nPerGroup:=.N,by=junctionID]
    dt[,idxGroup:=.GRP, by=junctionID]
    dt[,idxInGroup:=seq_len(.N),by=junctionID]
    
    if(max(dt$idxGroup) * m * freq < 10){
        freq <- 10/(max(dt$idxGroup) * m)
        warning("Injection-frequency is to low. Increasing it to `",
                signif(freq, 2), "` so we can inject at least 10 events ",
                "into the data set!")
    }
    
    # where do we inject and how
    indexOut_groups <- matrix(ncol=m, sample(c(0,1,-1),
            max(dt$idxGroup)*m, replace=TRUE, 
            prob=c(1-freq, freq/2, freq/2)))
        
    # positions where outliers will be injected
    list_index <- which(indexOut_groups != 0, arr.ind=TRUE)
    
    
    # sample primary injection
    primaryInjection <- merge(as.data.table(list_index), dt,
            sort=FALSE, by.x="row", by.y="idxGroup", allow.cartesian=TRUE)
    primaryInjection[,primary:=idxInGroup==sample(idxInGroup,1),by="row,col"]
    primaryInjection[,n:=n[cbind(idxInCount, col)]]
    primaryInjection[,k:=k[cbind(idxInCount, col)]]
    primaryInjection[,nOk:=length(unique(n)) == 1, by="row,col"]
    if(nrow(primaryInjection[nOk == FALSE]) / nrow(primaryInjection) > 0.05){
        # this can happen in the lower range if the prediction of the direction
        # did not work properly. Where the donor/acceptor goes the wrong way.
        # But this should not be too frequent. So ignore if it is below 5%
        warning("It looks like you have a problem with your indexes. ",
                "Please report it on github to the developer.")
    }
    primaryInjection <- primaryInjection[nOk == TRUE]
    primaryInjection[,I:=cumsum(primary)]
    primaryInjection[,I:=max(I),by="row,col"]
    
    
    # get matrix accessors (row, col) aka (junctionID, sampleID)
    primaryIndexInReal <- as.matrix(
        primaryInjection[primary == TRUE, .(idxInCount, col)])
    primaryIndexInGroup <- as.matrix(
        primaryInjection[primary == TRUE, .(row, col)])
    secondaryIndexInReal <- as.matrix(
        primaryInjection[primary != TRUE, .(idxInCount, col)])
    secondaryIndexInGroup <- as.matrix(
        primaryInjection[primary != TRUE, .(row, col)])
    # init with empty for psiSite. Function rep does not 
    # support times of length zero eg: rep(1, integer(0))
    secondaryPrimaryIndex <- as.matrix(primaryInjection[FALSE,.(
            row, col, idxInCount, junctionID, nPerGroup, I, rep=integer(.N))])
    if(nrow(primaryInjection[primary == TRUE & nPerGroup > 1]) > 0){
        secondaryPrimaryIndex <- as.matrix(primaryInjection[
                primary == TRUE & nPerGroup > 1, .(
                        idxInCount, junctionID, nPerGroup, I, 
                        rep=rep(TRUE, nPerGroup-1)),
                by="row,col"])
    }
    
    # get direction of injection
    injDirection <- indexOut_groups[primaryIndexInGroup]
    # table(injDirection)
    
    # get current psi of this junction and calculate maximal possible delta psi
    junction_psi <- psi[primaryIndexInReal]
    maxDpsi <- ifelse(injDirection > 0, 1 - junction_psi, junction_psi)
    
    # if not possible to inject here with at least minDpsi, 
    # change injection direction
    injectionPossible <- maxDpsi > minDpsi
    # table(injectionPossible)
    injDirection[!injectionPossible] <- -injDirection[!injectionPossible]
    maxDpsi <- ifelse(injDirection > 0, 1 - junction_psi, junction_psi)
    
    # ensure that injected points becomes an outlier by adding
    # mean delta psi as an offset
    meanDpsi <- rowMeans2(abs(psi[primaryIndexInReal[,"idxInCount"],] - 
            rowMeans2(psi[primaryIndexInReal[,"idxInCount"],])))
    
    # sample delta psi for injection from uniform distribution between 
    # min and max dpsi
    injMinDpsi <- ifelse(minDpsi + meanDpsi < maxDpsi,
                         minDpsi + meanDpsi, maxDpsi)
    injDpsi <- injDirection * switch(deltaDistr,
            uniformDistr = runif(length(injMinDpsi), injMinDpsi, maxDpsi),
            ifelse(deltaDistr > maxDpsi, maxDpsi, deltaDistr))
    
    # get N of this junction
    n_ji <- n[primaryIndexInReal]
    
    # inject new primary k_ij -> change k based on n and new_psi 
    # (and take pseudocounts into account): (k+1)/(n+2)=psi -> k = psi*(n+2)-1
    new_primary_psi <- junction_psi + injDpsi
    new_primary_k   <- pmax(pmin(round(
        new_primary_psi*(n_ji + 2*pseudocount()) - pseudocount()), n_ji), 0)
    
    # inject secondary k_ij -> change k based on n and the new_psi
    second_delta_psi <- - injDpsi[secondaryPrimaryIndex[,"I"]] * (
        psi[secondaryIndexInReal] / (
            1-psi[secondaryPrimaryIndex[,c("idxInCount", "col"),drop=FALSE]]))
    new_second_psi <- psi[secondaryIndexInReal] + second_delta_psi
    
    # sanity check for injected psi
    if(any(new_second_psi < -0.0001) | any(new_second_psi > 1.0001)){
        warning("Have to cut injected delta psi for: ", 
                sum(new_second_psi < -0.0001), " and ", 
                sum(new_second_psi > 1.0001), " instances!")
    }
    new_second_psi <- pmin(pmax(new_second_psi, 0), 1)
    new_second_k <- pmax(pmin(round(new_second_psi * ((
                    n[secondaryIndexInReal] + 2 * pseudocount()) - 
                            pseudocount())), 
            n[secondaryIndexInReal]), 0)
    
    # 
    # check if we really have an outlier or not
    # and get primiary/secondary injections
    # 
    goodInjections <- n[primaryIndexInReal] >= minCoverage &
        abs(k[primaryIndexInReal] - new_primary_k) > 2
    goodJunctionIds <- primaryInjection[primary == TRUE][
        goodInjections,junctionID]
    goodSecondary <- primaryInjection[primary == FALSE][,
            junctionID %in% goodJunctionIds]
    
    message(paste0(date(), ": Injecting outliers: ", sum(goodInjections), 
            " / ", sum(goodSecondary), " (primary/secondary)"))
    
    # 
    # prepare the injection into the dataset only for those which passed QC
    # 
    
    # set counts (k and o)
    k[primaryIndexInReal[goodInjections,,drop=FALSE]]  <- 
        new_primary_k[goodInjections]
    k[secondaryIndexInReal[goodSecondary,,drop=FALSE]] <- 
        new_second_k[goodSecondary]
    o <- n - k
    
    # set injection status (direction, primary, secondary)
    indexOut <- matrix(0, nrow=nrow(k), ncol=ncol(k))
    indexOut[primaryIndexInReal[goodInjections,]] <- 
            injDirection[goodInjections]
    indexOut[secondaryIndexInReal[goodSecondary,]] <- 
            2 * injDirection[secondaryPrimaryIndex[,"I"][goodSecondary]]
    
    # set injected delta psi
    indexDeltaPsi <- matrix(0, nrow=nrow(k), ncol=ncol(k))
    indexDeltaPsi[primaryIndexInReal[goodInjections,]] <- 
            injDpsi[goodInjections]
    indexDeltaPsi[secondaryIndexInReal[goodSecondary,]] <- 
            second_delta_psi[goodSecondary]
    
    # 
    # do the injection and save the additional informations in the object
    # 
    
    # copy original k and o
    if(type == "psiSite"){
        setAssayMatrix(fds, type="psiSite", "originalCounts") <-
            counts(fds, type="psiSite", side="ofInterest")
        setAssayMatrix(fds, type="psiSite", "originalOtherCounts") <-
            counts(fds, type="psiSite", side="other")
    } else {
        setAssayMatrix(fds, type=type, "originalCounts") <-
            counts(fds, type=type, side="ofInterest")
        setAssayMatrix(fds, type="psi5", "originalOtherCounts") <-
            counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "originalOtherCounts") <-
            counts(fds, type="psi3", side="other")
    }
    
    # store new k and o counts including the outlier counts
    counts(fds, type=type, side="ofInterest") <- k
    counts(fds, type=type, side="other") <- o
    
    # store positions of true outliers and their true delta PSI value
    setAssayMatrix(fds=fds, name="trueOutliers", type=type) <- indexOut
    setAssayMatrix(fds=fds, name="trueDeltaPSI", type=type) <- indexDeltaPsi
    
    return(fds)
}


removeInjectedOutliers <- function(fds, type){

    # copy injected k and o counts
    if(type == "psiSite"){
        setAssayMatrix(fds, type="psiSite", "outlierCounts") <- 
            counts(fds, type="psiSite", side="ofInterest")
        setAssayMatrix(fds, type="psiSite", "outlierOtherCounts") <- 
            counts(fds, type="psiSite", side="other")
    }
    else{
        setAssayMatrix(fds, type="psi5", "outlierCounts") <- 
            counts(fds, type="psi5", side="ofInterest")
        setAssayMatrix(fds, type="psi5", "outlierOtherCounts") <- 
            counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "outlierOtherCounts") <- 
            counts(fds, type="psi3", side="other")
    }

    # assign original k and o to rawCountsJ and rawOtherCounts
    if(type == "psiSite"){
        counts(fds, type="psiSite", side="ofInterest") <- 
            getAssayMatrix(fds, type="psiSite", "originalCounts")
        counts(fds, type="psiSite", side="other") <- 
            getAssayMatrix(fds, type="psiSite", "originalOtherCounts")

        assays(fds)[['originalCounts_psiSite']] <- NULL
        assays(fds)[['originalOtherCounts_psiSite']] <- NULL
    }
    else{
        counts(fds, type="psi5", side="ofInterest") <- 
            getAssayMatrix(fds, type="psi5", "originalCounts")
        counts(fds, type="psi5", side="other") <- 
            getAssayMatrix(fds, type="psi5", "originalOtherCounts")
        counts(fds, type="psi3", side="other") <-
            getAssayMatrix(fds, type="psi3", "originalOtherCounts")

        assays(fds)[['originalCounts_psi5']] <- NULL
        assays(fds)[['originalOtherCounts_psi5']] <- NULL
        assays(fds)[['originalOtherCounts_psi3']] <- NULL
    }

    return(fds)

}

restoreInjectedOutliers <- function(fds, type){

    # copy original k and o counts
    if(type == "psiSite"){
        setAssayMatrix(fds, type="psiSite", "originalCounts") <- 
            counts(fds, type="psiSite", side="ofInterest")
        setAssayMatrix(fds, type="psiSite", "originalOtherCounts") <- 
            counts(fds, type="psiSite", side="other")
    }
    else{
        setAssayMatrix(fds, type="psi5", "originalCounts") <- 
            counts(fds, type="psi5", side="ofInterest")
        setAssayMatrix(fds, type="psi5", "originalOtherCounts") <- 
            counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "originalOtherCounts") <- 
            counts(fds, type="psi3", side="other")
    }

    # assign injected k and o to rawCountsJ and rawOtherCounts
    if(type == "psiSite"){
        counts(fds, type="psiSite", side="ofInterest") <- 
            getAssayMatrix(fds, type="psiSite", "outlierCounts")
        counts(fds, type="psiSite", side="other") <- 
            getAssayMatrix(fds, type="psiSite", "outlierOtherCounts")

        assays(fds)[['outlierCounts_psiSite']] <- NULL
        assays(fds)[['outlierOtherCounts_psiSite']] <- NULL
    }
    else{
        counts(fds, type="psi5", side="ofInterest") <- 
            getAssayMatrix(fds, type="psi5", "outlierCounts")
        counts(fds, type="psi5", side="other") <- 
            getAssayMatrix(fds, type="psi5", "outlierOtherCounts")
        counts(fds, type="psi3", side="other") <-
            getAssayMatrix(fds, type="psi3", "outlierOtherCounts")

        assays(fds)[['outlierCounts_psi5']] <- NULL
        assays(fds)[['outlierOtherCounts_psi5']] <- NULL
        assays(fds)[['outlierOtherCounts_psi3']] <- NULL
    }

    return(fds)

}


# log2 fold change of measured psi vs AE predicted psi
log2fc <- function(realPsi, predictedPsi){
    return( log2(realPsi) - log2(predictedPsi) )
}
