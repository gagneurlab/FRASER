
# Plot methods
#
# plotVolcanoPerGene <- function(fds, sampleID, type=c("psi5", "psi3", "psiSite"),
#                             labelGenes=NULL, xoffset=0.8){
#     # extract data
#     dt2p <- data.table(
#         p=pVals(fds, type=type)[,sampleID],
#         z=getAssayMatrix(fds, "delta", type=type)[,sampleID],
#         geneID=mcols(fds, type=type)$hgnc_symbol)
#
#     # group by gene
#     dt2p2 <- dt2p[, p:=p.adjust(p, method="holm"), by=geneID]
#     dt2p2 <- dt2p2[order(geneID, p)]  [!duplicated(geneID) & !is.na(geneID)]
#     dt2p2 <- dt2p2[, padj:=p.adjust(p, method="BY")]
#     maxPByFDR <- dt2p2[padj < 0.1, max(p)]
#
#     # get x label
#     xlab <- switch(type,
#                    'psi3' = bquote(Delta * Psi[3]),
#                    'psi5' = bquote(Delta * Psi[5]),
#                    'psiSite' = bquote(Delta ~ "SE")
#     )
#
#     # plot it
#     g <- ggplot(dt2p2, aes(x=z, y=-log10(p))) + geom_point(col="gray70", alpha=0.4) +
#         xlab(xlab) +
#         ylab(bquote(-log[10]~"(p value)")) +
#         geom_vline(xintercept=c(-0.3, 0.3), color="firebrick", linetype=2) +
#         geom_hline(yintercept=-log10(maxPByFDR), color="firebrick", linetype=4) +
#         ggtitle(paste("Volcano plot:", sampleID, " - ", type)) +
#         geom_point(data=dt2p2[abs(z) > 0.3 & padj < 0.1], aes(x=z, y=-log10(p)), col="firebrick")
#
#     for(lg in labelGenes){
#         g <- g + dt2p2[geneID == lg,
#                        annotate("text", x=z*xoffset, y=-log10(p)*1.0, label=geneID)]
#     }
#     g
# }

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
                          x=expression(-log[10] ~  "(exp. p value)"),
                          y=expression(-log[10] ~ "(obs. p value)") ) +
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
                          x=expression(-log[10] ~  "(exp. p value)"),
                          y=expression(-log[10] ~ "(obs. p value)") ) +
        geom_ribbon(alpha=0.2, col="gray", aes(x=exp, ymin = lower, ymax = upper)) +
        geom_segment(aes(x=min(exp), y=min(exp), xend=max(exp), yend=max(exp)),
                     col="firebrick")
    g

}


plotExpression <- function(fds, type=c("psi5", "psi3", "psiSite"),
                    site=NULL, result=NULL, colGroup=NULL, basePlot=TRUE,
                    title=paste0("Expression plot: ", site), ...){
    if(!is.null(result)){
        type <- as.character(result$type)
        site <- getIndexFromResultTable(fds, result)
    } else {
        type <- match.arg(type)
    }

    dt <- getPlottingDT(fds, axis="row", type=type, idx=site, ...)

    if(!is.null(colGroup)){
        if(all(colGroup %in% samples(fds))){
            colGroup <- samples(fds) %in% colGroup
        }
        dt[colGroup,aberrant:=TRUE]
    }

    # # estimate point densities
    # # adapted from http://auguga.blogspot.com/2015/10/r-heat-scatter-plot.html
    # dens <- kde2d(dt$ni,dt$ki,
    #               h = c(ifelse(bandwidth.nrd(dt$ni) == 0, 0.1, bandwidth.nrd(dt$ni)),
    #                     ifelse(bandwidth.nrd(dt$ki) == 0, 0.1, bandwidth.nrd(dt$ki))))
    # # create a new data frame of that 2d density grid
    # gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
    # names(gr) <- c("xgr", "ygr", "zgr")
    # # Fit a model
    # mod <- loess(zgr~xgr*ygr, data=gr)
    # # Apply the model to the original data to estimate density at that point
    # dt[,pointdens:=predict(mod, newdata=data.frame(xgr=dt$ni, ygr=dt$ki))]
    # g <- ggplot(dt, aes(x=ni, y=ki)) + geom_point(aes(color=pointdens),
    #                                               show.legend=FALSE) +
    #     scale_color_gradientn(colors = colorpalette('heat', 5)) +

    g <- ggplot(dt, aes(x=n, y=k, color=!aberrant, text=paste0(
                "Sample: ", sampleID, "<br>",
                "Counts (K): ", k, "<br>",
                "Total counts (N): ", n, "<br>",
                "p value: ", signif(pval, 5), "<br>",
                "padjust: ", signif(padj, 5), "<br>",
                "Observed Psi: ", round(obsPsi, 2), "<br>",
                "Predicted mu: ", round(predPsi), "<br>"))) +
        geom_point(alpha=ifelse(dt$aberrant, 1, 0.7)) +
        scale_x_log10() +
        scale_y_log10() +
        geom_abline(intercept=0, slope=1) +
        theme_bw() +
        theme(legend.position="none") +
        xlab("Total junction coverage + 2 (N)") +
        ylab("Junction count + 1 (K)") +
        ggtitle(title)

    plotBasePlot(g, basePlot)
}

plotExpectedVsObservedPsi <- function(fds, type=c("psi5", "psi3", "psiSite"),
                    idx=NULL, result=NULL, colGroup=NULL, main=NULL,
                    basePlot=TRUE, padjCutoff=0.05){
    type <- match.arg(type)

    # get plotting data
    dt <- getPlottingDT(fds, axis="row", result=result, idx=idx)
    idx <- unique(dt$idx)

    if(is.null(main)){
        # TODO extract feature name if present: siteName <- ...
        main <- paste0("Expression vs prediction plot: ", idx)
    }

    dt[,aberrant:=FALSE]
    dt[padj < padjCutoff, aberrant:=TRUE]

    if(!is.null(colGroup)){
        if(is.logical(colGroup)){
            dt[colGroup, aberrant:=TRUE]
        } else {
            warning("not implemented yet!")
        }
    }

    xlab <- switch(type,
        'psi3' = bquote("Observed " ~ psi[3]),
        'psi5' = bquote("Observed " ~ psi[5]),
        'psiSite' = "Observed SE"
    )
    ylab <- switch(type,
        'psi3' = bquote("Predicted " ~ psi[3]),
        'psi5' = bquote("Predicted " ~ psi[5]),
        'psiSite' = "Predicted SE"
    )

    g <- ggplot(dt, aes(x=obsPsi, y=predPsi, color=!aberrant)) +
        geom_point(alpha=ifelse(dt$aberrant, 1, 0.5)) +
        geom_abline(intercept = 0, slope=1) +
        theme_bw() +
        theme(legend.position="none") +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(main)

    plotBasePlot(g, basePlot)
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
        scale_y_log10() +
        labs(title=main, x="Sample rank", y="Number of outliers")

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
        scale_y_log10() +
        labs(title=main, x="Sample rank", y="Number of outliers") +
        scale_color_discrete(name="", labels = c(bquote(Psi[3]), bquote(Psi[5]), "SE"))

    g

}

plotQQPerJunction <- function(fds, type=c("psi5", "psi3", "psiSite"), site,
                          maxOutlier=2, conf.alpha=0.05, breakTies=TRUE,
                          sample=FALSE, highlightOutliers=FALSE, cutYaxis=FALSE,
                          mainName=paste0("QQ-Plot (", type, ")")){

    # extract data
    pvals <- pVals(fds, type, byGroup=FALSE)[site,]
    padj  <- padjVals(fds, type, byGroup=FALSE)[site,]

    outlier <- logical(length(pvals))
    if(isTRUE(highlightOutliers)){
        outlier[padj < 0.1] <- TRUE
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
    g <- ggplot(dt[plotPoint,], aes(x=exp, y=obs, color=outlier)) +
        geom_point(show.legend=FALSE) +
        scale_color_manual(values=c("black", "firebrick")) +
        theme_bw() + labs(title=mainName,
                          x=expression(-log[10] ~  "(exp. p value)"),
                          y=expression(-log[10] ~ "(obs. p value)") ) +
        geom_ribbon(alpha=0.2, col="gray", aes(x=exp, ymin = lower, ymax = upper)) +
        geom_segment(aes(x=min(exp), y=min(exp), xend=max(exp), yend=max(exp)),
                     col="firebrick")
    g

}
