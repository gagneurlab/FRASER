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

    # store information about the simulation in the fds
    # (same values for psi3 and psi5)
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
    junctionGroups <- sort(sample(x = donorGroups, size=j, replace = TRUE))   # assign junction to donor/acceptor groups

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
    # Simulate n for each junction using lognormal and negative binomial
    #
    junction_n <- rlnorm(j, meanlog=nbMeanLog, sdlog=nbSize)
    n <- t(vapply(junction_n, function(x){
        # draw counts for every sample from negative binomial
        n <- rnbinom(m, mu=x, size=nbSize)},
        double(m)))
    # n <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)

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
    group_rho <- rlnorm(d, meanlog=rhoMeanlog, sdlog=rhoSdlog)  # group dispersion
    # group_rho <- abs(rnorm(d, mean=0.0001, sd=0.05))     # group dispersion
    dt1 <- data.table(group=donorGroups, rho=group_rho)
    dt2 <- data.table(group=junctionGroups)
    rhoFull <- merge(dt1, dt2, by="group")
    rho <- rhoFull[,rho]                           # set dispersion for every junction



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
        psi_mu <- group_mu[-1,]
        if(is.null(nrow(psi_mu))){ # only one junction in the group, convert to matrix with one row so that colSums works
            psi_mu <- matrix(psi_mu, nrow=1)
        }
        # scale mu's so that they sum up to 1 again (after mu for SE is removed)
        psi_mu    <- psi_mu /  matrix(colSums(psi_mu), nrow=nrow(psi_mu), ncol=ncol(psi_mu), byrow = TRUE)

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
#' @rdname makeSimulatedFraserDataSet
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
            simulatedPSI = getAssayMatrix(fds, "truePSI", type=type) )
    
    # 
    # matrix of donor/acceptor for possible injection
    # 
    if(type == "psiSite"){
        dt <- as.data.table(rowRanges(fds, type))[,.(
                chr=factor(seqnames), start, end, strand=factor(strand),
                junctionID=spliceSiteID)]
    }else{
        getInjectionSite <- function(type, strand, startId, endId){
            selectionMat <- as.matrix(data.frame(row=seq_along(startId), col=
                    1 + as.vector(
                            type == "psi5" & strand == "-" | 
                            type == "psi3" & strand == "+")))
            ans <- as.matrix(cbind(startId, endId, strand))
            ans[selectionMat]
        }
        dt <- as.data.table(rowRanges(fds, type))
        dt <- dt[,.(chr=factor(seqnames), start, end, strand=factor(strand),
                junctionID=getInjectionSite(type, strand, startID, endID))]
    }
    
    # Get groups where outlier can be injected
    available_groups <- dt[,.(.I, junctionID)][,
            .(nPerGroup=.N, idxGroup=.GRP, idxInGroup=seq_len(.N), 
                    idxInCount=I),by=junctionID][order(idxInCount)]
    
    if(type %in% c("psi5", "psi3")){
        if(nrow(available_groups[nPerGroup > 1]) > 0){
            available_groups <- available_groups[nPerGroup > 1]
            available_groups[,idxGroup:=.GRP,by=junctionID]
        } else {
            warning("No alternative splicing detected in the dataset. ",
                    "Use full dataset to inject arteficial outliers.")
        }
    }
    
    if(max(available_groups$idxGroup) * m * freq < 10){
        freq <- 10/(max(available_groups$idxGroup) * m)
        warning("Injection-frequency is to low. Increasing it to `",
                signif(freq, 2), "` so we can inject at least 10 events ",
                "into the data set!")
    }
    
    # where do we inject and how
    indexOut_groups <- matrix(ncol=m, sample(c(0,1,-1),
            max(available_groups$idxGroup)*m, replace=TRUE, 
            prob=c(1-freq, freq/2, freq/2)))
        
    # positions where outliers will be injected
    list_index <- which(indexOut_groups != 0, arr.ind=TRUE)
    
    
    # sample primary injection
    primaryInjection <- merge(as.data.table(list_index), available_groups,
            sort=FALSE, by.x="row", by.y="idxGroup", allow.cartesian=TRUE)
    primaryInjection[,primary:=idxInGroup==sample(idxInGroup,1),by="row,col"]
    primaryInjection[,n:=n[cbind(idxInCount, col)]]
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
            1-psi[secondaryPrimaryIndex[,c("idxInCount", "col")]]))
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
            " / ", sum(goodSecondary), " (primary/secondary"))
    
    # 
    # prepare the injection into the dataset only for those which passed QC
    # 
    
    # set counts (k and o)
    k[primaryIndexInReal[goodInjections,]]  <- new_primary_k[goodInjections]
    k[secondaryIndexInReal[goodSecondary,]] <- new_second_k[goodSecondary]
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
        second_delta_psi[secondaryPrimaryIndex[,"I"][goodSecondary]]
    
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

debug_old_code <- function(){
    # apply injection function to each outlier
    message(date(), ": Injecting ", nrow(list_index), " outliers ...")
    result <- bplapply(seq_row(list_index), list_index=list_index,
            available_groups=available_groups, deltaDistr=deltaDistr,
            indexOut_groups=indexOut_groups, type=type, psi=psi, n=n, dt=dt,
            minDpsi=minDpsi, verbose=verbose, BPPARAM=BPPARAM,
            FUN=injectSingleOutlier)
    
    
    # copy original k and o
    if(type == "psiSite"){
        setAssayMatrix(fds, type="psiSite", "originalCounts") <-
            counts(fds, type="psiSite", side="ofInterest")
        setAssayMatrix(fds, type="psiSite", "originalOtherCounts") <-
            counts(fds, type="psiSite", side="other")
    }
    else{
        setAssayMatrix(fds, type=type, "originalCounts") <-
            counts(fds, type=type, side="ofInterest")
        setAssayMatrix(fds, type="psi5", "originalOtherCounts") <-
            counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "originalOtherCounts") <-
            counts(fds, type="psi3", side="other")
    }
    
    # get all junctions, samples, ... where outlier injection changed the counts
    junctions     <- unlist(lapply(result, "[[", 'junctions'))
    samples       <- unlist(lapply(result, "[[", 'sample'))
    newKs         <- unlist(lapply(result, "[[", 'newKs'))
    newOs         <- unlist(lapply(result, "[[", 'newOs'))
    injDirection  <- unlist(lapply(result, "[[", 'injDirections'))
    injDeltaPSI   <- unlist(lapply(result, "[[", 'injDeltaPSI'))

    # matrices to store indices of outliers and their delta psi
    indexOut <- matrix(0, nrow = j, ncol = m)
    indexDeltaPSI <- matrix(0, nrow = j, ncol = m)

    # set counts to changed counts after the injection
    replaceIndices                <- matrix(c(junctions,samples), ncol=2)
    k[replaceIndices]             <- newKs
    if(length(available_groups) == j){
        o                           <- n - k
    } else{
        o[replaceIndices]           <- newOs
    }
    indexOut[replaceIndices]      <- injDirection
    indexDeltaPSI[replaceIndices] <- injDeltaPSI

    # store positions of true outliers and their true delta PSI value
    setAssayMatrix(fds=fds, name="trueOutliers", type=type) <- indexOut
    setAssayMatrix(fds=fds, name="trueDeltaPSI", type=type) <- indexDeltaPSI

    # store new k counts which include the outlier counts
    counts(fds, type=type, side="ofInterest") <- k

    # store modified other counts
    counts(fds, type=type, side="other") <- o

    return(fds)

}

injectSingleOutlier <- function(j, available_groups, list_index, deltaDistr, 
                    indexOut_groups, type, psi, n, dt=dt, minDpsi, verbose){
    
    # extract group, sample and injecetion direction (i.e +1/up or -1/down)
    row       <- list_index[j,'row']
    group     <- available_groups[row]
    sample    <- list_index[j,'col']
    injDirection    <- indexOut_groups[row, sample]
    
    # sample one junction from all junction within this group
    group_junctions <- which(dt[, junctionID == group])
    junction <- sample(group_junctions, 1)
    
    # get current psi of this junction and calculate maximla possible value of delta psi for the injection
    junction_psi    <- psi[junction, sample]
    maxDpsi <- junction_psi
    if(injDirection > 0){
        maxDpsi <- 1 - junction_psi
    }
    
    # if not possible to inject here with at least minDpsi, change injection direction
    if(maxDpsi < minDpsi){
        injDirection <- -injDirection
        indexOut_groups[row, sample] <- injDirection
        maxDpsi <- junction_psi
        if(injDirection > 0){
            1 - junction_psi
        }
    }
    
    # ensure that injected points becomes an outlier by adding mean delta psi as an offset
    meanDpsi <- mean(abs(psi[junction,] - mean(psi[junction,])))
    
    # sample delta psi for injection from uniform distribution between min and max dpsi
    minDpsi <- ifelse(minDpsi+meanDpsi < maxDpsi, minDpsi+meanDpsi, maxDpsi)
    injDpsi <- injDirection * switch(deltaDistr,
            uniformDistr = runif(1, minDpsi, maxDpsi),
            ifelse(as.double(deltaDistr) > maxDpsi, maxDpsi, as.double(deltaDistr)) )
    
    # if(is.na(injDpsi)){
    #   browser()
    # }
    
    # get N of this junction
    n_j <- n[group_junctions,sample]
    
    # new counts after injection
    newKs     <- integer(length(group_junctions))
    indexDpsi <- double(length(group_junctions))
    indOut    <- integer(length(group_junctions))
    
    # for all other junctions in this group
    for(i in seq_along(group_junctions)){
        junction_k <- group_junctions[i]
        # get new_psi and change k based on n and new_psi (and take pseudocounts into account): (k+1)/(n+2)=psi -> k = psi*(n+2) - 1
        if(junction_k == junction){
            new_psi <- junction_psi + injDpsi
            new_k <- round( new_psi*(n_j[i] + 2*pseudocount()) - pseudocount() )
            
            # store position of outlier
            indOut[i]      <- injDirection
            indexDpsi[i]   <- injDpsi
        }
        else{
            deltaPSI_k <- - (psi[junction_k,sample] / (1-junction_psi)) * injDpsi
            new_psi <- min(1, max(0, psi[junction_k,sample] + deltaPSI_k))
            new_k <- round( new_psi*(n_j[i] + 2*pseudocount()) - pseudocount() )
            
            # also store position of other junction used in swap
            indOut[i]      <- -injDirection * 2
            indexDpsi[i]   <- deltaPSI_k
        }
        
        # ensure that k is smaller than n and not negative
        newKs[i] <- max(0, min(new_k, n_j[i]))
    }
    
    if(verbose){
        print(paste("injected outlier", j, "with delta PSI of", 
                injDpsi, "at junction", junction, "and sample", sample))
    }
    
    return(list(newKs = newKs, newOs=sum(newKs)-newKs, 
            junctions = group_junctions, injDeltaPSI = indexDpsi, 
            injDirections = indOut, 
            sample = rep(sample, length(group_junctions))))
}


#
# Inject artificial outliers in an existing fds
#
injectOutliersBySwapping <- function(fds, type=type, nrOutliers=500, deltaPSI=pmin(0.7, pmax(-0.7, rnorm(nrOutliers, 0, 0.5))),
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

        possible_junctions <- if(deltaPSI[i] > 0){ which(rowMeans(psi) - 0.2 < (1 - deltaPSI[i])) }else{ which(rowMeans(psi) + 0.2 > -deltaPSI[i]) }
        # possible_junctions <- ifelse(deltaPSI[i] > 0, which(rowMaxs(psi) < (1 - deltaPSI[i])), which(rowMaxs(psi) > -deltaPSI[i]) )
        counter <- 0

        # for each outlier, draw random junction-sample pair and check if outlier can be injected there
        while(!successful && counter <= 1000){

            counter <- counter + 1
            # if(verbose){
            #   print(paste0(counter, ". try"))
            # }

            if(length(intersect(available_junctions, possible_junctions)) == 0 ){
                break
            }

            junction <- sample(intersect(available_junctions, possible_junctions), 1)
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
                    dt <- data.table(
                        junctionID = 1:j,
                        chr = as.factor(seqnames(fds)),
                        start = start(fds),
                        end = end(fds),
                        strand = as.factor(strand(fds))
                    )
                    if(type == "psi5"){
                        samePos <- dt[start == dt[junctionID == junction, start] & chr == dt[junctionID == junction, chr] & strand == dt[junctionID == junction, strand], junctionID]
                    }
                    else if(type == "psi3"){
                        samePos <- dt[end == dt[junctionID == junction, end] & chr == dt[junctionID == junction, chr] & strand == dt[junctionID == junction, strand], junctionID]
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
                            # indexDeltaPSI[junction, sample] <- deltaPSI[i]
                            indexDeltaPSI[group_others[g], sample] <- ( (k[group_others[g],sample] + pseudocount())/(n[group_others[g],sample] + 2*pseudocount()) ) - psi[group_others[g], sample] # new dPsi - old dPsi

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

        if(successful){

            # remove junction where outlier was injected from available junctions (each junction at most one outlier)
            available_junctions <- available_junctions[available_junctions != junction]

            if(verbose){
                print(paste("Injected outlier", i, "at index", junction, ",", sample, ": new k =", k_new, ", n =", n_ji))
            }

        }
        else{

            if(verbose){
                print(paste("Skipped outlier", i, "with delta psi =", deltaPSI[i]))
            }

        }


    }

    # store positions of true outliers and their true delta PSI value
    setAssayMatrix(fds=fds, name="trueOutliers", type=type) <- indexOut
    setAssayMatrix(fds=fds, name="trueDeltaPSI", type=type) <- indexDeltaPSI

    # store new k counts which include the outlier counts
    counts(fds, type=type, side="ofInterest") <- k

    # for SE: also store modified junction counts
    if(type == "psiSite"){
        counts(fds, type="psi5", side="ofInterest") <- k_other
    }

    # re-calculate other counts if swapping was done
    if(swap){
        # if psi values were already calculated, delete them (or rather move the old values somewhere else?) so that they are recalculated
        for(psiType in c("psi3", "psi5", "psiSite")){
            if(assayExists(fds, psiType)){
                assays(fds)[[psiType]] <- NULL
            }
            if(assayExists(fds, paste0('delta_', psiType))){
                assays(fds)[[paste0('delta_', psiType)]] <- NULL
            }
        }
        fds <- calculatePSIValues(fds, overwriteCts = TRUE)
    }
    else{
        # when no swapping was done, just put other counts as n-k
        counts(fds, type=type, side="other") <- n-k
    }

    return(fds)

}

removeInjectedOutliers <- function(fds, type){

    # copy injected k and o counts
    if(type == "psiSite"){
        setAssayMatrix(fds, type="psiSite", "outlierCounts") <- counts(fds, type="psiSite", side="ofInterest")
        setAssayMatrix(fds, type="psiSite", "outlierOtherCounts") <- counts(fds, type="psiSite", side="other")
    }
    else{
        setAssayMatrix(fds, type="psi5", "outlierCounts") <- counts(fds, type="psi5", side="ofInterest")
        setAssayMatrix(fds, type="psi5", "outlierOtherCounts") <- counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "outlierOtherCounts") <- counts(fds, type="psi3", side="other")
    }

    # assign original k and o to rawCountsJ and rawOtherCounts
    if(type == "psiSite"){
        counts(fds, type="psiSite", side="ofInterest") <- getAssayMatrix(fds, type="psiSite", "originalCounts")
        counts(fds, type="psiSite", side="other") <- getAssayMatrix(fds, type="psiSite", "originalOtherCounts")

        assays(fds)[['originalCounts_psiSite']] <- NULL
        assays(fds)[['originalOtherCounts_psiSite']] <- NULL
    }
    else{
        counts(fds, type="psi5", side="ofInterest") <- getAssayMatrix(fds, type="psi5", "originalCounts")
        counts(fds, type="psi5", side="other") <- getAssayMatrix(fds, type="psi5", "originalOtherCounts")
        counts(fds, type="psi3", side="other") <-getAssayMatrix(fds, type="psi3", "originalOtherCounts")

        assays(fds)[['originalCounts_psi5']] <- NULL
        assays(fds)[['originalOtherCounts_psi5']] <- NULL
        assays(fds)[['originalOtherCounts_psi3']] <- NULL
    }

    return(fds)

}

restoreInjectedOutliers <- function(fds, type){

    # copy original k and o counts
    if(type == "psiSite"){
        setAssayMatrix(fds, type="psiSite", "originalCounts") <- counts(fds, type="psiSite", side="ofInterest")
        setAssayMatrix(fds, type="psiSite", "originalOtherCounts") <- counts(fds, type="psiSite", side="other")
    }
    else{
        setAssayMatrix(fds, type="psi5", "originalCounts") <- counts(fds, type="psi5", side="ofInterest")
        setAssayMatrix(fds, type="psi5", "originalOtherCounts") <- counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "originalOtherCounts") <- counts(fds, type="psi3", side="other")
    }

    # assign injected k and o to rawCountsJ and rawOtherCounts
    if(type == "psiSite"){
        counts(fds, type="psiSite", side="ofInterest") <- getAssayMatrix(fds, type="psiSite", "outlierCounts")
        counts(fds, type="psiSite", side="other") <- getAssayMatrix(fds, type="psiSite", "outlierOtherCounts")

        assays(fds)[['outlierCounts_psiSite']] <- NULL
        assays(fds)[['outlierOtherCounts_psiSite']] <- NULL
    }
    else{
        counts(fds, type="psi5", side="ofInterest") <- getAssayMatrix(fds, type="psi5", "outlierCounts")
        counts(fds, type="psi5", side="other") <- getAssayMatrix(fds, type="psi5", "outlierOtherCounts")
        counts(fds, type="psi3", side="other") <-getAssayMatrix(fds, type="psi3", "outlierOtherCounts")

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
