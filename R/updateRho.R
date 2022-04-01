#'
#' Update theta step for autoencoder
#'
#' @noRd
updateRho <- function(fds, type, rhoRange, BPPARAM, verbose){
    
    k <- K(fds)
    n <- N(fds)
    y <- predictY(fds, noiseAlpha=currentNoiseAlpha(fds))
    
    # fitparameters <- bplapply(seq_len(nrow(k)), estRho, nll=truncNLL_rho,
    #         k=k, n=n, y=y, rhoRange=rhoRange, BPPARAM=BPPARAM)
    fitparameters <- bplapply(seq_len(nrow(k)), estRho, nll=trunc_negLogLikelihoodRho_penalized,
                              k=k, n=n, y=y, rhoRange=rhoRange, lambda=0,
                              BPPARAM=BPPARAM)
    
    rho(fds) <- plogis(vapply(fitparameters, "[[", 
            double(1), "minimum"))
    
    if(isTRUE(verbose)){
        stxt <- capture.output(summary(rho(fds)))
        message(date(), ": rho fit:\n\t", paste(stxt, collapse="\n\t"))
    }
    
    validObject(fds)
    return(fds)
}

estRho <- function(idx, k, n, y, rhoRange, nll, control=list(), lambda=0){
    ki <- k[idx,]
    ni <- n[idx,]
    yi <- y[idx,]
    
    # est <- optimize(f=nll, interval=rhoRange, yi=yi, ki=ki, ni=ni, 
    #                 maximum=FALSE, tol=0.0000001)
    # est
    est <- optimize(f=nll, interval=c(-40, 40), mui=plogis(yi), ki=ki, ni=ni, lambda=lambda,
                    maximum=FALSE, tol=0.0000001)
    est
}

negLogLikelihoodRho <- function(rho, ki, ni, mui){
    #-mean(dbetabinom(ki + 0.5, ni + 1, mu, rho, log=TRUE))
    
    r  <- (1-rho)/rho
    eps <- 0.5
    alpha  <- lgamma(mui*r)
    alphaK <- lgamma(mui*r + ki + eps)
    beta   <- lgamma((mui-1)*(-r))
    betaNK <- lgamma((mui-1)*(-r) + (ni - ki + eps))
    
    #mean negative log likelihood with pseudocounts
    mean(alpha + beta - alphaK - betaNK - lgamma(ni+1+2*eps) + 
            lgamma(ki+1+eps) + lgamma(ni-ki+1+eps) + lgamma(r + ni + 2*eps) - 
            lgamma(r))
}

trunc_negLogLikelihoodRho <- function(rho, ki, ni, mui){
    #-mean(dbetabinom(ki + 0.5, ni + 1, mu, rho, log=TRUE))
    
    r  <- (1-rho)/rho
    eps <- 0.5
    alpha  <- lgamma(mui*r)
    alphaK <- lgamma(mui*r + ki + eps)
    beta   <- lgamma((mui-1)*(-r))
    betaNK <- lgamma((mui-1)*(-r) + (ni - ki + eps))
    
    #mean negative log likelihood with pseudocounts
    mean(alpha + beta - alphaK - betaNK )
}

trunc_negLogLikelihoodRho_penalized <- function(logit_rho, ki, ni, mui, lambda){
    #-mean(dbetabinom(ki, ni, mui, rho, log=TRUE))
    
    rho <- plogis(logit_rho)
    r  <- (1-rho)/rho
    alpha  <- lgamma(mui*r)
    alphaK <- lgamma(mui*r + ki)
    beta   <- lgamma((mui-1)*(-r))
    betaNK <- lgamma((mui-1)*(-r) + (ni - ki))
    
    #mean negative log likelihood with pseudocounts
    mean(alpha + beta - alphaK - betaNK ) + lambda * (logit_rho*logit_rho)
}


methodOfMomentsRho <- function(k, n, rhoRange=c(1e-5, 1 - 1e-5)){
    # taken from wiki: 
    # https://en.wikipedia.org/wiki/Beta-binomial_distribution#Method_of_moments
    k <- as.matrix(k) + pseudocount()
    mode(k) <- "double"
    n <- as.matrix(n) + 2*pseudocount()
    mode(n) <- "double"
    
    # normalizing counts with respect to mean(n) because for each junction an 
    # equal n for all samples is needed
    nM <- round(rowMeans(n))
    kNorm <- round(nM/n * k)
    
    N  <- ncol(kNorm)
    m1 <- rowSums(kNorm)/N
    m2 <- rowSums(kNorm*kNorm)/N
    
    denom <- (nM*(m2/m1 - m1 - 1) + m1)
    a <- (nM*m1 - m2) / denom
    b <- ((nM - m1)*(nM - m2/m1))/ denom
    rho <- 1 / (1 + a + b)
    
    # set correct boundaries
    rho[a < 0 | b < 0] <- rhoRange[1]
    rho[is.na(rho)] <- rhoRange[1]
    rho <- pmax(rho, rhoRange[1])
    rho <- pmin(rho, rhoRange[2])
    
    return(rho)
}
