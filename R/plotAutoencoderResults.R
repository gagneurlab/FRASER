
#plot junctions after autoencoder fit
plotJunction <- function(idx, fds, type=currentType(fds), threshold=5e-2, dist="BetaBinomial", paper=FALSE){
    par.old <- par(no.readonly=TRUE)
    on.exit(par(par.old))

    if(isTRUE(paper)){
        par(mfrow=c(2,2))
    } else {
        par(mfrow=c(2,3))
    }

    K <- K(fds)
    N <- N(fds)
    mu <- predictedMeans(fds)
    padj  <- padjVals(fds, dist=dist)

    plotData(idx, K, N, padj, threshold)
    #plotPsi(idx, fds, K, N, padj, threshold)
    #plotMu(idx, mu, fds, padj, threshold)
    plotPsiVsMu(idx, K, N, mu, padj, threshold)

    if(!is.null(pVals(fds, dist=dist))){
        if(isFALSE(paper)){
            plotPvalHist(idx, fds, dist)
        }
        plotQQ(idx, fds, "BetaBinomial")
        if(isFALSE(paper)){
            plotQQ(idx, fds, "Binomial")
        }
    }
    plotZscore(idx, fds)
}

plotData <- function(idx, K, N, pvals, threshold=5e-2, main=NULL){
    if(is.null(main)){
        main <- paste("Idx:", idx)
    }
    heatscatter(N[idx,] + 2*pseudocount(), K[idx,] + pseudocount(),
            xlab="N (Junction Coverage)", ylab="K (Reads of interest)", pch=16,
            main=main, log="xy", cexplot = 1)

    if(!missing(pvals) & !is.null(pvals)){
        pval = pvals[idx,]
        pos <- which(pval < threshold)
        points(N[idx,pos] + 2*pseudocount(), K[idx, pos] + pseudocount(), pch=19, col="firebrick")
    }
    grid()
    abline(0,1)
}

plotPsiVsMu <- function(idx, K, N, mu, pvals, threshold=5e-2){
    psi <- (K + pseudocount())/(N + (2*pseudocount()))
    plot(y=psi[idx,], x=mu[idx,], xlab="predictedMu", ylab="raw PSI", pch=19)
    grid()
    abline(0, 1, col="firebrick")
    if(!missing(pvals) & !is.null(pvals)){
        pval = pvals[idx,]
        pos <- which(pval < threshold)
        points(y=psi[idx,pos], x=mu[idx,pos], pch=2, col="green")
    }
}

plotPsi <- function(idx, fds, K, N, pvals, threshold=5e-2){
    psi <- (K[idx,] + pseudocount())/(N[idx,] + (2*pseudocount()))
    plot(psi, xlab="sample", ylab="raw PSI", pch=19)

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

plotPvalHist <- function(idx, fds, dist){
    hist(pVals(fds, dist=dist)[idx,], xlab="P-values", main=paste("Pvalues (", dist, ") for", idx))
}

plotQQ <- function(idx, fds, dist, main=NULL, threshold=5e-2){
    pval <- pVals(fds, dist=dist)[idx,]
    exp <- -log10(ppoints(length(pval)))
    obs <- -log10(sort(pval))
    if(is.null(main)){
        main <- paste("QQ-plot", "(", dist, ") for", idx)
    }

    xlim=range(exp)
    ylim=range(obs)

    plot(NA, xlim=xlim, ylim=ylim,
       main=main,
       xlab=expression(
           paste(-log[10], " (expected ", italic(P), "-value)")),
       ylab=expression(
           paste(-log[10], " (observed ", italic(P), "-value)")))

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
    points(exp, obs, pch=16)

    if(any(obs > -log10(threshold))){
        pos <- which(obs > -log10(threshold))
        points(exp[pos], obs[pos], pch=19, col="firebrick")
    }

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



plotLoss <- function(fds, type){
    type
    fds <- fds_new
    lossList <- metadata(fds)[[paste0('loss_', type)]]
    dt2plot <- as.data.table(melt(do.call(rbind, list(
        max=colMaxs(lossList, na.rm=TRUE),
        mean=colMeans(lossList, na.rm=TRUE),
        min=colMins(lossList, na.rm=TRUE),
        sd=colSds(lossList, na.rm=TRUE))), value.name="Loss"))
    colnames(dt2plot)[1:2] <- c("Aggregation", "IterSteps")

    # set iterations
    dt2plot[grepl("init",IterSteps),iteration:=0]
    dt2plot[is.na(iteration),iteration:=as.numeric(as.character(
        gsub("_.*", "" ,gsub("final_", "" , IterSteps))))]

    # set step
    dt2plot[grepl("final", IterSteps),Step:="D fit"]
    dt2plot[is.na(Step) & !grepl("init", IterSteps),Step:="E fit"]
    dt2plot[is.na(Step),Step:="Init"]
    dt2plot[,Iteration:=1:.N/1,by="Aggregation,Step"]
    dt2plot[Step=="E fit", Iteration:=Iteration/3]
    dt2plot[Step=="D fit", Iteration:=Iteration/2]

    # summaries
    ggplot(dt2plot[Aggregation != "sd"], aes(x=Iteration, y=Loss, col=Step, lty=Aggregation)) +
        geom_ribbon(data=dt2plot[Aggregation == "mean",], aes(
            ymin = dt2plot[Aggregation == "mean", Loss] - dt2plot[Aggregation == "sd", Loss],
            ymax = dt2plot[Aggregation == "mean", Loss] + dt2plot[Aggregation == "sd", Loss]), fill="gray80", col=NA) +
        geom_line() +
        scale_y_log10()
}
