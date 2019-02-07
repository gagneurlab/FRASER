
#plot junctions after autoencoder fit
plotJunction <- function(idx, fds, threshold=5e-2, dist="BetaBinom"){
    par.old <- par(no.readonly=TRUE)
    on.exit(par(par.old))

    par(mfrow=c(2,3))

    K <- K(fds)
    N <- N(fds)
    mu <- predictedMeans(fds)
    pvals <- pVals(fds, dist=dist)
    padj  <- padjVals(fds, dist=dist)

    plotData(idx, K, N, padj, threshold)
    plotPsi(idx, fds, K, N, padj, threshold)
    plotMu(idx, mu, fds, padj, threshold)

    if(!is.null(pvals)){
        plotPvalHist(idx, pvals)
        plotQQ(idx, pvals, "BetaBinomial")
        plotQQ(idx, fds, "Binomial")
    }
    #plotZscore(idx, fds)
}

plotData <- function(idx, K, N, pvals, threshold=5e-2){
    heatscatter(N[idx,] + (2*pseudocount()), K[idx,] + pseudocount(),
            xlab="N (Junction Coverage)", ylab="K (Reads of interest)", pch=16,
            main=paste("Idx:", idx), log="xy", cexplot = 1)

    if(!missing(pvals) & !is.null(pvals)){
        pval = pvals[idx,]
        pos <- which(pval < threshold)
        points(N[idx,pos] + 1, K[idx, pos] + 0.5, pch=2, col="green")
    }
    grid()
    abline(0,1)
}

plotPsi <- function(idx, fds, K, N, pvals, threshold=5e-2){
    psi <- (K[idx,] + pseudocount())/(N[idx,] + (2*pseudocount()))
    plot(1:ncol(fds), psi, xlab="sample", ylab="raw PSI", pch=19)

    if(!missing(pvals) & !is.null(pvals)){
        pval = pvals[idx,]
        pos <- which(pval < threshold)
        points(pos, psi[pos], pch=2, col="green")
    }
}

plotMu <- function(idx, mu, fds, pvals, threshold=5e-2){

    plot(1:ncol(fds), sort(mu[idx,]), xlab="sample", ylab="predicted mu",
            main = paste("rho:", signif(rho(fds)[idx], 5)), pch=19)

    if(!missing(pvals) & !is.null(pvals)){
        pval = pvals[idx,]
        pos <- rank(mu)[which(pval < threshold)]
        points(pos, mu[idx,pos], pch=2, col="green")
    }
}

plotPvalHist <- function(idx, pvals){
    hist(pvals[idx,], xlab="P-values", main=paste("Pvalues for", idx))
}

plotQQ <- function(idx, pvals, distribution){
    pval <- pvals[idx,]
    exp <- -log10(ppoints(length(pval)))
    obs <- -log10(sort(pval))
    main <- paste("QQ-plot", "(",distribution,") for", idx)

    xlim=range(exp)
    ylim=range(obs)

    plot(NA, xlim=xlim, ylim=ylim,
       main=main, xlab="expected p-val", ylab="observered")

    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    conf.alpha <- 0.05

    if(is.numeric(conf.alpha)){
    len <- length(exp)
    slen <- seq_len(len)
    getY <- function(x, exp){
      x1 <- exp[2]
      x2 <- exp[1]
      y1 <- -log10(x[2])
      y2 <- -log10(x[1])
      m <- (y2-y1)/(x2-x1)
      return(10^-(y1 + m*((x2+1)-x1)))
    }
    upper <- qbeta(    conf.alpha/2, slen, rev(slen))
    lower <- qbeta(1 - conf.alpha/2, slen, rev(slen))
    polygon(col="gray", border="gray", x=c(rev(exp), max(exp)+c(1,1), exp),
            y=-log10(c(
              rev(upper), getY(upper, exp), getY(lower, exp), lower)))
    }

    # Add points
    points(exp, obs)

    # diagonal and grid
    abline(0,1,col="firebrick")
    grid()

}

# plotQQbinomial <- function(idx, pvals){#fds, mu){
#   # pval <- singlePvalueBinomial(idx, K(fds), N(fds), mu)
#   pval <- pvals[idx,]
#   main <- paste("QQ-plot (binomial) for", idx)
#
#   plot(-log10(ppoints(length(pval))), -log10(sort(pval)),
#        main=main, xlab="expected p-val", ylab="observered")
#   grid()
#   abline(0,1)
# }

plotZscore <- function(idx, fds){
  z <- zScores(fds)[idx,]
  hist(z)
}


