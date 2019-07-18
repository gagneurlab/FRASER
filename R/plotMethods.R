
# Plot methods

plotVolcanoPlot <- function(fds, sampleID, type=c("psi5", "psi3", "psiSite"),
                            labelGenes=NULL){
    # extract data
    dt2p <- data.table(
        p=pVals(fds, type=type)[,sampleID],
        padj=padjVals(fds, type=type)[,sampleID],
        z=getAssayMatrix(fds, "delta", type=type)[,sampleID],
        geneID=mcols(fds, type=type)$hgnc_symbol)

    # group by gene
    dt2p2 <- dt2p[, p:=p.adjust(p, method="holm"), by=geneID]
    dt2p2 <- dt2p2[order(geneID, p)]  [!duplicated(geneID) & !is.na(geneID)]
    maxPByFDR <- dt2p2[padj < 0.1, max(p)]

    # plot it
    g <- ggplot(dt2p2, aes(x=z, y=-log10(p))) + geom_point() +
        xlab("delta PSI") +
        ylab("-log10(P-value)") +
        geom_vline(xintercept=c(-0.3, 0.3), color="firebrick", linetype=2) +
        geom_hline(yintercept=-log10(maxPByFDR), color="firebrick", linetype=4) +
        ggtitle(paste("Volcano plot:", sampleID, " - ", type))

    for(lg in labelGenes){
        g <- g + dt2p2[geneID == lg,
                       annotate("text", x=z, y=-log10(p)*1.0, label=geneID)]
    }
    g
}

plotQQPlot <- function(fds, type=c("psi5", "psi3", "psiSite"), site, sampleID=NULL,
                       maxOutlier=2, conf.alpha=0.05, breakTies=TRUE,
                       sample=FALSE, highlightOutliers=TRUE,
                       mainName=paste0("QQ-Plot (", type, ")"), cutYaxis=FALSE){

    # extract data
    pvals <- getPvalsPerGene(fds, type, sampleID=sampleID)

    # points
    obs <- -log10(sort(pvals))
    obs[is.na(obs)] <- 0
    if(length(unique(obs)) < 2 | length(obs) < 2){
        warning("There are no pvalues or all are NA!")
        return(FALSE)
    }
    if(breakTies){
        obs <- breakTies(obs, logBase=10, decreasing=TRUE)
    }
    exp <- -log10(ppoints(length(obs)))
    len <- length(exp)

    # limits for nice plotting
    maxPoint <- max(c(exp, obs))
    ylim <- range(0, min(exp[1]*maxOutlier, maxPoint))
    outlier <- FALSE
    if(isTRUE(highlightOutliers) & !is.null(sampleID)){
        outlier <- obs > max(ylim)
    }


    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(is.numeric(conf.alpha)){
        getY <- function(x, exp){
            x1 <- exp[2]
            x2 <- exp[1]
            y1 <- -log10(x[2])
            y2 <- -log10(x[1])
            m <- (y2-y1)/(x2-x1)
            return(10^-(y1 + m*((x2+1)-x1)))
        }
        upper <- qbeta(conf.alpha/2,   1:len, rev(1:len))
        lower <- qbeta(1-conf.alpha/2, 1:len, rev(1:len))
    }

    plotPoint <- TRUE
    if(isTRUE(sample)){
        lo <- length(obs)
        plotPoint <- 1:lo %in% unique(c(1:min(lo, 100), sort(sample(1:lo,
                                                                    size = min(lo, 30000), prob=log10(1+lo:1)/sum(log10(1+lo:1))))))
    }

    # get plot data
    dt <- data.table( obs=obs, exp=exp, lower=-log10(lower), upper=-log10(upper),
                      plotPoint=plotPoint, outlier=outlier)
    if(isTRUE(cutYaxis)){
        dt[outlier,obs:=rep(max(ylim), sum(outlier))]
    }

    # create qq-plot
    g <- ggplot(dt[plotPoint,], aes(x=exp, y=obs)) + geom_point(aes(color=outlier),
                                                                show.legend=FALSE) +
        scale_color_manual(values=c(1,2)) +
        theme_bw() + labs(title=mainName,
                          x=expression(-log[10] ~  "(expected P-value)"),
                          y=expression(-log[10] ~ "(observed P-value)") ) +
        geom_ribbon(alpha=0.2, col="gray", aes(x=exp, ymin = lower, ymax = upper)) +
        geom_segment(aes(x=min(exp), y=min(exp), xend=max(exp), yend=max(exp)),
                     col="firebrick")
    g

}


plotJunctionCounts <- function(fds, type=c("psi5", "psi3", "psiSite"), site){

    k <- K(fds, type)
    n <- N(fds, type)

    dt <- data.table(ki=k[site,], ni=n[site,])

    g <- ggplot(dt, aes(x=(ni+2), y=(ki+1))) + geom_point() +
        scale_x_log10() + scale_y_log10() +
        geom_abline(intercept = 0, slope=1) + theme_bw() +
        xlab("N (Total Site Coverage) + 2") + ylab("K (Junction Counts) + 1")
    g

}

plotOutlierSampleRank <- function(fds, type=c("psi5", "psi3", "psiSite"),
                                  main=paste0("Number outliers per sample (", type, ")"),
                                  alpha=0.05){

    pvals <- getPvalsPerGene(fds, type, pvals=padjVals(fds, type))

    nrOutliers <- apply(pvals, 2, function(x){sum(x < alpha)})

    dt <- data.table(nrHits=nrOutliers)
    setorder(dt, nrHits)
    dt[,rank:=1:.N]

    g <- ggplot(dt, aes(x=rank, y=nrHits)) + geom_line() + theme_bw() +
        labs(title=main, x="Sample Rank", y="Number of Outliers")

    g

}
