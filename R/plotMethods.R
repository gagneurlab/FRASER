
# Plot methods

plotVolcanoPerGene <- function(fds, sampleID, type=c("psi5", "psi3", "psiSite"),
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

plotQQPerGene <- function(fds, type=c("psi5", "psi3", "psiSite"), gene=NULL,
                          sampleID=NULL, maxOutlier=2, conf.alpha=0.05,
                          breakTies=TRUE, sample=FALSE, highlightOutliers=TRUE,
                          mainName=paste0("QQ-Plot (", type, ")"), cutYaxis=FALSE){

    # extract data
    pvals <- getPvalsPerGene(fds, type, sampleID=sampleID)
    if(!is.null(gene)){
        index <- which(getGeneIDs(fds, type) == gene)
        pvals <- pvals[index,]
    }

    # points
    obs <- -log10(sort(c(pvals)))
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
    if(isTRUE(highlightOutliers) & (!is.null(gene) & !is.null(sampleID))){
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
        scale_color_manual(values=c("black", "firebrick")) +
        theme_bw() + labs(title=mainName,
                          x=expression(-log[10] ~  "(expected P-value)"),
                          y=expression(-log[10] ~ "(observed P-value)") ) +
        geom_ribbon(alpha=0.2, col="gray", aes(x=exp, ymin = lower, ymax = upper)) +
        geom_segment(aes(x=min(exp), y=min(exp), xend=max(exp), yend=max(exp)),
                     col="firebrick")
    g

}


plotGlobalQQPerGene <- function(fds, types=psiTypes, maxOutlier=2, conf.alpha=0.05,
                          breakTies=TRUE, sample=TRUE, mainName="QQ-Plot",
                          BPPARAM=MulticoreParam(min(bpworkers(), length(types))) ){

    dtls <- bplapply(types, function(type){
        # extract data
        pvals <- getPvalsPerGene(fds, type, sampleID=NULL)

        # points
        obs <- -log10(sort(c(pvals)))
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
            plotPoint <- 1:lo %in% unique(c(1:min(lo, 1000), sort(sample(1:lo,
                                                                         size = min(lo, 30000), prob=log10(1+lo:1)/sum(log10(1+lo:1))))))
        }

        # get plot data
        dt <- data.table( obs=obs, exp=exp, lower=-log10(lower), upper=-log10(upper),
                          plotPoint=plotPoint, type=type)
        dt
    }, BPPARAM=BPPARAM)

    dt2p <- rbindlist(dtls)
    dt2p <- dt2p[plotPoint == TRUE,]

    # create qq-plot
    g <- ggplot(dt2p, aes(x=exp, y=obs, color=type)) + geom_point() +
        scale_color_discrete(name="", labels=c(bquote(Psi[3]), bquote(Psi[5]), "SE")) +
        labs(title=mainName,
                          x=expression(-log[10] ~  "(expected P-value)"),
                          y=expression(-log[10] ~ "(observed P-value)") ) +
        geom_ribbon(alpha=0.2, col="gray", aes(x=exp, ymin = lower, ymax = upper)) +
        geom_segment(aes(x=min(exp), y=min(exp), xend=max(exp), yend=max(exp)),
                     col="firebrick")
    g

}


plotJunctionCounts <- function(fds, type=c("psi5", "psi3", "psiSite"), site,
                               highlightSample=NULL,
                               title=paste0("Site: ", site)){

    k <- K(fds, type)
    n <- N(fds, type)

    dt <- data.table(ki=(k[site,]+pseudocount()), ni=(n[site,]+2*pseudocount()))

    # estimate point densities
    # adapted from http://auguga.blogspot.com/2015/10/r-heat-scatter-plot.html
    dens <- kde2d(dt$ni,dt$ki,
                  h = c(ifelse(bandwidth.nrd(dt$ni) == 0, 0.1, bandwidth.nrd(dt$ni)),
                        ifelse(bandwidth.nrd(dt$ki) == 0, 0.1, bandwidth.nrd(dt$ki))))
    # create a new data frame of that 2d density grid
    gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    # Apply the model to the original data to estimate density at that point
    dt[,pointdens:=predict(mod, newdata=data.frame(xgr=dt$ni, ygr=dt$ki))]

    g <- ggplot(dt, aes(x=ni, y=ki)) + geom_point(aes(color=pointdens),
                                                  show.legend=FALSE) +
        scale_x_log10() + scale_y_log10() +
        scale_color_gradientn(colors = colorpalette('heat', 5)) +
        geom_abline(intercept = 0, slope=1) + theme_bw() +
        xlab("N (Total Junction Coverage) + 2") + ylab("K (Junction Counts) + 1") +
        ggtitle(title)
    if(!is.null(highlightSample)){
        g <- g + geom_point(data=dt[highlightSample, ], aes(x=ni, y=ki), colour="green", shape=17, size=3)
    }
    g

}

plotOutlierSampleRankPerGenePerType <- function(fds, type=c("psi5", "psi3", "psiSite"),
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

plotOutlierSampleRankPerGene <- function(fds, alpha=0.05, types=psiTypes, main=paste0(
                                                "Number outliers per sample"),
                                                BPPARAM=MulticoreParam(
                                                    min(bpworkers(), length(types))) ){

    dtls <- bplapply(types, function(type){
        pvals <- getPvalsPerGene(fds, type, pvals=padjVals(fds, type))
        nrOutliers <- apply(pvals, 2, function(x){sum(x < alpha)})
        dt <- data.table(type=type, nrHits=nrOutliers)
        setorder(dt, nrHits)
        dt[,rank:=1:.N]
        dt
    }, BPPARAM=BPPARAM)

    dt2p <- rbindlist(dtls)

    g <- ggplot(dt2p, aes(x=rank, y=nrHits, color=type)) + geom_line() +
        labs(title=main, x="Sample Rank", y="Number of Outliers") +
        scale_color_discrete(name="", labels = c(bquote(Psi[3]), bquote(Psi[5]), "SE"))

    g

}

plotQQPerJunction <- function(fds, type=c("psi5", "psi3", "psiSite"), site,
                          maxOutlier=2, conf.alpha=0.05, breakTies=TRUE,
                          sample=FALSE, highlightOutliers=NULL, cutYaxis=FALSE,
                          mainName=paste0("QQ-Plot (", type, ")")){

    # extract data
    pvals <- pVals(fds, type, byGroup=FALSE)[site,]

    outlier <- rep(FALSE, length(c(pvals)))
    if(!is.null(highlightOutliers)){
        outlier[highlightOutliers] <- TRUE
    }

    # points
    sortOrder <- order(c(pvals))
    outlier <- outlier[sortOrder]
    obs <- -log10(sort(c(pvals)))
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
        scale_color_manual(values=c("black", "firebrick")) +
        theme_bw() + labs(title=mainName,
                          x=expression(-log[10] ~  "(expected P-value)"),
                          y=expression(-log[10] ~ "(observed P-value)") ) +
        geom_ribbon(alpha=0.2, col="gray", aes(x=exp, ymin = lower, ymax = upper)) +
        geom_segment(aes(x=min(exp), y=min(exp), xend=max(exp), yend=max(exp)),
                     col="firebrick")
    g

}
