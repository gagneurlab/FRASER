
#plot junctions after autoencoder fit
plotJunction <- function(idx, fds, pvals=pVals(fds), mu=predictMu(fds), threshold=1e-2, setpar=TRUE){
  if(isTRUE(setpar)){
    par(mfrow=c(2,3))
  }
  
  K <- K(fds)
  N <- N(fds)
  
  padj <- NULL
  if(!is.null(pvals)){
    padj <- p.adjust(pvals, "BH")
    padj <- matrix(padj, nrow=dim(fds)[1], ncol=dim(fds)[2])
  }
  
  plotData(idx, K, N, padj, threshold)
  plotPsi(idx, fds, K, N, padj, threshold)
  plotMu(idx, mu, fds, padj, threshold)
  
  if(!is.null(pvals)){
    plotPvalHist(idx, pvals)
    plotQQ(idx, pvals)
  }
  plotZscore(idx, fds)
  
}

plotData <- function(idx, K, N, pvals=NULL, threshold=1e-2){
  heatscatter(N[idx,] + 1, K[idx,] + 0.5, xlab="N (Junction Coverage)", ylab="K (Reads of interest)", pch=16,
              main=paste("Idx:", idx), log="xy", cexplot = 1)
  if(!is.null(pvals)){
    pval = pvals[idx,]
    pos <- which(pval < threshold)
    points(N[idx,pos] + 1, K[idx, pos] + 0.5, pch=2, col="green")
  }
  grid()
  abline(0,1)
}

plotPsi <- function(idx, fds, K, N, pvals, threshold=1e-2){
  psi <- (K[idx,] + 0.5)/(N[idx,] + 1)
  plot(1:ncol(fds), psi, xlab="sample", ylab="raw PSI", pch=19)
  if(!is.null(pvals)){
    pval = pvals[idx,]
    pos <- which(pval < threshold)
    points(pos, psi[pos], pch=2, col="green")
  }
}

plotMu <- function(idx, mu, fds, pvals, threshold=1e-2){
  
  plot(1:ncol(fds), mu[idx,], xlab="sample", ylab="predicted mu", main = paste("rho:", signif(rho(fds)[idx], 5)), pch=19)
  if(!is.null(pvals)){
    pval = pvals[idx,]
    pos <- which(pval < threshold)
    points(pos, mu[idx,pos], pch=2, col="green")
  }
  
}

plotPvalHist <- function(idx, pvals){
  hist(pvals[idx,], xlab="P-values", main=paste("Pvalues for", idx))
}

plotQQ <- function(idx, pvals){
  pval <- pvals[idx,]
  main <- paste("QQ-plot for", idx)
  
  plot(-log10(ppoints(length(pval))), -log10(sort(pval)),
       main=main, xlab="expected p-val", ylab="observered")
  grid()
  abline(0,1)
}


plotZscore <- function(idx, fds){
  z <- zScores(fds)[idx,]
  hist(z)
}


