#' Visualization functions for FraseR
#'
#' The FraseR package provides mutliple functions to visualize
#' the data and the results of a full data set analysis.
#'
#' This is the list of all plotting function provided by FraseR:
#' \itemize{
#'   \item plotAberrantPerSample()
#'   \item plotVolcano()
#'   \item plotExpression()
#'   \item plotQQ()
#'   \item plotExpectedVsObservedPsi()
#'   \item plotCountCorHeatmap()
#'   \item plotFilterExpression()
#'   \item plotEncDimSearch()
#' }
#'
#' For a detailed description of each plot function please see the details.
#' Most of the functions share the same parameters.
#'
#### Data specific parameters
#' @param fds An FraseRDataSet object.
#' @param type The psi type: either psi5, psi3 or psiSite (for SE).
#' @param sampleID A sample ID which should be plotted. Can also be a vector.
#'             Integers are treated as indices.
#' @param idx,site A junction site ID or gene ID or one of both, which
#'             should be plotted. Can also be a vector. Integers are treated
#'             as indices.
#' @param padjCutoff,zScoreCutoff,deltaPsiCutoff Significance, Z-score or delta
#'             psi cutoff to mark outliers
#' @param global Flag to plot a global Q-Q plot, default FALSE
#' @param normalized If TRUE, the normalized psi values are used, the default,
#'             otherwise the raw psi values
#' @param aggregate,aggregated If TRUE, the pvalues are aggregated by gene, 
#'             otherwise junction level pvalues are used (default).
#' @param result The result table to be used by the method.
#' @param BPPARAM BiocParallel parameter to use.
#' @param Ncpus Number of cores to use.
#### Graphical parameters
#' @param main Title for the plot, if missing a default title will be used.
#' @param colGroup Group of samples that should be colored.
#' @param basePlot if \code{TRUE} (default), use the R base plot version, else
#'             use the plotly framework.
#' @param conf.alpha If set, a confidence interval is plotted, defaults to 0.05
#' @param samplePoints Sample points for Q-Q plot, defaults to max 30k points
#' @param logit If TRUE, the default, psi values are plotted in logit space.
#' @param nClust Number of clusters to show in the row and
#'             column dendrograms.
#' @param sampleClustering A clustering of the samples that should be used as an
#'             annotation of the heatmap.
#' @param show_rownames,show_colnames Logical value indicating whether to show
#'             row or column names on the heatmap axes.
#' @param annotation_col,annotation_row Row or column annotations that should be
#'             plotted on the heatmap.
#' @param plotType Character string indicating the type of correlation heatmap
#'             that should be plotted. Can be either 'sampleCorrelation' for a
#'             sample-sample correlation heatmap or 'junctionSample' for a
#'             junction-sample correlation heatmap.
#' @param topN,topJ Top x most variable junctions that should be used in the
#'             heatmap. TopN is used for sample-sample correlation heatmaps and
#'             topJ for junction-sample correlation heatmaps.
#' @param minMedian,minDeltaPsi Minimal median value or minimal delta psi of a
#'             junction to be considered for the correlation heatmap.
#' @param border_color Sets the border color of the heatmap
#' @param plotMeanPsi,plotCov If \code{TRUE}, then the heatmap is annotated with
#'             the mean psi values or the junction coverage.
#' @param bins Set the number of bins to be used in the histogram.
#' @param legend.position Set legend position (x and y coordinate), defaults to
#'             the top right corner.
#'
#### Additional ... parameter
#' @param ... Additional parameters passed to plot() or plot_ly() if not stated
#'             otherwise in the details for each plot function
#'
#' @details
#'
#' \code{plotAberrantPerSample}: The number of aberrant events per sample are
#' plotted sorted by rank. The ... parameters are passed on to the
#' \code{\link{aberrant}} function.
#'
#' \code{plotVolcano}: the volcano plot is sample-centric. It plots for a given
#' sample and psi type the negative log10 nominal P-values against the delta psi
#' values for all splice sites or aggregates by gene if requested.
#'
#' \code{plotExpression}: This function plots for a given site the
#' read count at this site (i.e. K) against the total coverage (i.e. N) for the
#' given psi type (psi5, psi3 or SE (psiSite)) for all samples.
#'
#' \code{plotQQ}: the quantile-quantile plot for a given gene or if
#' \code{global} is set to \code{TRUE} over the full data set. Here the
#' observed P-values are plotted against the expected ones in the negative
#' log10 space.
#'
#' \code{plotExpectedVsObservedPsi}: A scatter plot of the observed psi
#' against the predicted psi for a given site.
#'
#' \code{plotCountCorHeatmap}: The correlation heatmap of the count data either
#' of the full data set (i.e. sample-sample correlations) or of the top x most
#' variable junctions (i.e. junction-sample correlations). By default the values
#' are log transformed and row centered. The ... arguments are passed to the
#' \code{\link[pheatmap]{pheatmap}} function.
#'
#' \code{plotFilterExpression}: The distribution of FPKM values. If the 
#' FraseRDataSet object contains the \code{passedFilter} column, it will plot 
#' both FPKM distributions for the expressed genes and for the filtered genes.
#'
#' \code{plotEncDimSearch}: Visualization of the hyperparameter optimization.
#' It plots the encoding dimension against the achieved loss (area under the
#' precision-recall curve). From this plot the optimum should be choosen for
#' the \code{q} in fitting process.
#'
#' @return If base R graphics are used nothing is returned else the plotly or
#'             the gplot object is returned.
#'
#' @name plotFunctions
#' @rdname plotFunctions
#' @aliases plotFunctions plotAberrantPerSample plotVolcano plotQQ 
#'             plotExpression plotCountCorHeatmap plotFilterExpression 
#'             plotExpectedVsObservedPsi plotEncDimSearch
#' @examples
#' # TODO
#' TODO <- 1
#'
NULL

#'
#' Volcano plot
#'
#' Plots the p values over the delta psi values, known as volcano plot.
#' Visualizes per sample the outliers. By type and aggregate by
#' gene if requested.
#'
#' @rdname plotFunctions
#' @export
plotVolcano <- function(fds, sampleID, type, deltaPsiCutoff=0.3,
                    padjCutoff=0.05, basePlot=TRUE, aggregate=FALSE,
                    main=paste0("Volcano plot: ", sampleID)){

    dt <- getPlottingDT(fds, axis="col", type=type, idx=sampleID,
            aggregate=aggregate)

    padj_line <- dt[padj < 0.05, min(5.5, -log10(pval))]
    if(length(padj_line) == 0){
        padj_line <- 5.5
    }

    dt[,aberrant:=FALSE]
    dt[abs(deltaPsi) >= deltaPsiCutoff & padj <= padjCutoff, aberrant:=TRUE]

    g <- ggplot(dt, aes(x=deltaPsi, y=-log10(pval), color=aberrant, text=paste0(
        "SampleID: ", sampleID, "<br>",
        "featureID: ", featureID, "<br>",
        "Raw count (K): ", k, "<br>",
        "Raw total count (N): ", n, "<br>",
        "p value: ", signif(pval, 5), "<br>",
        "delta Psi: ", round(deltaPsi, 2), "<br>",
        "Type: ", type))) +
        geom_point(aes(alpha=ifelse(isTRUE(aberrant), 1, 0.5))) +
        xlab(expression(paste(Delta, Psi))) +
        ylab(expression(paste(-log[10], "(p value)"))) +
        geom_vline(xintercept=c(-deltaPsiCutoff, deltaPsiCutoff),
                   color="firebrick", linetype=2) +
        geom_hline(yintercept=padj_line, color="firebrick", linetype=4) +
        ggtitle(main) +
        theme_bw() +
        theme(legend.position="none") +
        scale_color_manual(values=c("gray70", "firebrick"))

    if(isFALSE(basePlot)){
        g <- g + xlab("delta Psi") +
            ylab("-log[10](p value)")
    }
    plotBasePlot(g, basePlot)
}

#'
#' Number of aberrant events per sample
#'
#' Plot the number of aberrant events per samples
#'
#' @rdname plotFunctions
#' @export
plotAberrantPerSample <- function(fds, main, type=c("psi3", "psi5", "psiSite"),
                    padjCutoff=0.1, zScoreCutoff=NA, deltaPsiCutoff=0.3,
                    aggregated=TRUE, BPPARAM=bpparam(), ...){

    type <- match.arg(type, several.ok=TRUE)

    if(missing(main)){
        main <- paste('Aberrant events per sample')
        if(isTRUE(aggregated)){
            main <- paste(main, "by gene")
        }
    }

    # extract outliers
    outliers <- bplapply(type, aberrant, fds=fds, by="sample",
            padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
            deltaPsiCutoff=deltaPsiCutoff, ..., BPPARAM=BPPARAM)
    dt2p <- rbindlist(lapply(seq_along(outliers), function(idx){
        vals <- outliers[[idx]]
        data.table(type=type[idx], value=sort(vals), median=median(vals),
                rank=seq_along(vals))
    }))

    # plot them
    g <- ggplot(dt2p, aes(x=rank, y=value, color=type)) +
        geom_line() +
        geom_hline(aes(yintercept=median, color=type, lty="Median")) +
        theme_bw() +
        theme_cowplot() +
        ggtitle(main) +
        xlab("Sample rank") +
        ylab("Number of outliers") +
        scale_color_brewer(palette="Dark2", name="Splice metric",
                labels=ggplotLabelPsi(dt2p[,unique(type)])) +
        scale_linetype_manual(name="", values=2, labels="Median")

    if(!all(dt2p[,value] == 0)){
        g <- g + scale_y_log10()
    }
    
    g
}


#'
#' Junction expression plot
#'
#' Plots the observed split reads of the junction of interest over all reads
#' coming from the given donor/acceptor.
#'
#' @rdname plotFunctions
#' @export
plotExpression <- function(fds, type=c("psi5", "psi3", "psiSite"),
                           site=NULL, result=NULL, colGroup=NULL, basePlot=TRUE,
                           main=NULL, ...){
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
    dt[,aberrant:=factor(aberrant, levels=c("TRUE", "FALSE"))]

    if(is.null(main)){
        main <- as.expression(bquote(bold(paste(
                .(ggplotLabelPsi(type)[[1]]), " expression plot: ",
                bolditalic(.(as.character(dt[,unique(featureID)]))),
                " (site ", .(as.character(dt[,unique(idx)])), ")"))))
    }

    g <- ggplot(dt, aes(x=n + 2, y=k + 1, color=aberrant, text=paste0(
            "Sample: ", sampleID, "<br>",
            "Counts (K): ", k, "<br>",
            "Total counts (N): ", n, "<br>",
            "p value: ", signif(pval, 5), "<br>",
            "padjust: ", signif(padj, 5), "<br>",
            "Observed Psi: ", round(obsPsi, 2), "<br>",
            "Predicted mu: ", round(predPsi), "<br>"))) +
        geom_point(alpha=ifelse(as.character(dt$aberrant) == "TRUE", 1, 0.7)) +
        scale_x_log10() +
        scale_y_log10() +
        theme_bw() +
        theme(legend.position="none", title=) +
        xlab("Total junction coverage + 2 (N)") +
        ylab("Junction count + 1 (K)") +
        ggtitle(main)

    if(is.null(colGroup)){
        g <- g + scale_colour_manual(
                values=c("FALSE"="gray70", "TRUE"="firebrick"))
    }

    plotBasePlot(g, basePlot)
}


#'
#' Expected over Overserved plot
#'
#' Plots the expected psi value over the observed psi value of the given
#' junction.
#'
#' @rdname plotFunctions
#' @export
plotExpectedVsObservedPsi <- function(fds, type=c("psi5", "psi3", "psiSite"),
                                      idx=NULL, result=NULL, colGroup=NULL,
                                      main=NULL, basePlot=TRUE, 
                                      padjCutoff=0.05){
    type <- match.arg(type)

    # get plotting data
    dt <- getPlottingDT(fds, axis="row", type=type, result=result, idx=idx)
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


#'
#' Q-Q plot
#'
#' Plots the quantile-quantile plot
#'
#' @rdname plotFunctions
#' @export
plotQQ <- function(fds, type=NULL, idx=NULL, result=NULL, aggregate=FALSE,
                    global=FALSE, main=NULL, conf.alpha=0.05,
                    samplePoints=30000, basePlot=TRUE,
                    Ncpus=min(4, getDTthreads()), ...){

    # check parameters
    if(is.null(aggregate)){
        aggregate <- isTRUE(global)
    } else if(!(is.logical(aggregate) |
                all(aggregate %in% colnames(mcols(fds))))){
        stop("Please provide TRUE/FALSE or a ",
             "charactor matching a column name in mcols.")
    }

    if(isTRUE(global)){
        if(is.null(type)){
            type <- psiTypes
        }
        dt <- rbindlist(lapply(type, getPlottingDT, fds=fds, axis="col",
                idx=TRUE, aggregate=aggregate, mc.cores=Ncpus, ...))
    } else {
        dt <- getPlottingDT(fds, axis="row", type=type, idx=idx, result=result,
                aggregate=aggregate, ...)
    }
    if(is.null(main)){
        if(isTRUE(global)){
            main <- "Global QQ plot"
        } else {
            type <- as.character(dt[,unique(type)])
            featureID <- as.character(dt[,unique(featureID)])
            main <- as.expression(bquote(bold(paste(
                    .(ggplotLabelPsi(type)[[1]]),
                    " Q-Q plot: ", bolditalic(.(featureID)),
                    " (site ", .(as.character(dt[,unique(idx)])), ")"))))
        }
    }


    # points
    dt2p <- dt[order(type, pval)]
    dt2p[,obs:=-log10(pval)]
    dt2p[is.na(obs), obs:=0]
    dt2p[is.infinite(obs), obs:=dt2p[is.finite(obs),max(obs)]]
    if(dt2p[,length(unique(obs))] < 2 | nrow(dt2p) < 2){
        warning("There are no pvalues or all are NA or 1!")
        return(FALSE)
    }
    dt2p[,exp:=-log10(ppoints(.N)),by=type]

    # down sample if requested
    dt2p[,plotPoint:=TRUE]
    if(isTRUE(samplePoints) | isScalarNumeric(samplePoints)){
        if(isScalarLogical(samplePoints)){
            samplePoints <- 30000
        }
        getProb <- function(n){
            rnl <- log10(rev(seq_len(n)) + 1)
            rnl/sum(rnl)
        }
        dt2p[,plotPoint:=
                rep(c(TRUE, FALSE), c(min(10000, .N), .N - min(10000, .N))) |
                sample(c(TRUE, FALSE), size=.N, replace=TRUE,
                        prob=c(samplePoints/.N, 1)), by=type]
    }

    # create qq-plot
    g <- ggplot(dt2p[plotPoint == TRUE,], aes(x=exp, y=obs, col=aberrant,
            text=paste(
                "<br>SampleID: ", sampleID, "<br>K: ", k, "<br>N: ", n))) +
        geom_point() +
        theme_bw() +
        theme(legend.position="none") +
        ggtitle(main) +
        xlab(expression(-log[10]~"(exp. p value)")) +
        ylab(expression(-log[10]~"(obs. p value)"))

    # Set color scale for global/local
    if(isFALSE(global)){
        g <- g + scale_color_manual(values=c("black", "firebrick"),
                name="Aberrant")
    } else {
        g$mapping$colour <- quote(type)
        g <- g + scale_color_brewer(palette="Dark2", name="Splice metric",
                labels=ggplotLabelPsi(unique(dt2p$type)))
    }


    # add confidence band if requesded
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(is.numeric(conf.alpha)){
        dt2p[,upper:=-log10(
                qbeta(conf.alpha/2,   seq_len(.N), rev(seq_len(.N)))), by=type]
        dt2p[,lower:=-log10(
                qbeta(1-conf.alpha/2, seq_len(.N), rev(seq_len(.N)))), by=type]
        g <- g + geom_ribbon(data=dt2p[plotPoint == TRUE,],
                aes(x=exp, ymin=lower, ymax=upper, text=NULL),
                alpha=0.2, color="gray")
    }

    # add abline in the end
    g <- g + geom_abline(intercept=0, slope=1, col="firebrick")
    if(isFALSE(global)){
        return(plotBasePlot(g, basePlot))
    }
    g
}


#'
#' Plots the results from the hyperparamter optimization.
#'
#' @rdname plotFunctions
#' @export
plotEncDimSearch <- function(fds, type=c("psi3", "psi5", "psiSite")){
    type <- match.arg(type)
    data <- hyperParams(fds, type=type, all=TRUE)
    if (is.null(data)) {
        warning(paste("no hyperparameters were estimated for", type, 
                      "\nPlease use `optimHyperParams` to compute them."))
        return()
    }
    if(!"nsubset" %in% colnames(data)){
        data[,nsubset:="NA"]
    }
    data[,noise:=as.factor(noise)]
    data[,nsubset:=as.factor(nsubset)]

    g1 <- ggplot(data, aes(q, aroc, col=nsubset, linetype=noise)) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        ggtitle("Q estimation") +
        xlab("Estimated q") +
        ylab("Area under the ROC curve")

    g2 <- ggplot(data, aes(q, eval, col=nsubset, linetype=noise)) +
        geom_point() +
        geom_line() +
        geom_smooth() +
        ggtitle("Q estimation") +
        xlab("Estimated q") +
        ylab("Model loss")


    ggarrange(g1, g2, nrow=2)
}


#'
#' Plot filter expression
#'
#' Histogram of the geometric mean per junction based on the filter status
#'
#' @rdname plotFunctions
#' @export
plotFilterExpression <- function(fds, bins=200, legend.position=c(0.8, 0.8)){
    cts    <- K(fds, "psi5")
    rowlgm <- exp(rowMeans(log(cts + 1)))

    dt <- data.table(
            value=rowlgm,
            passed=mcols(fds, type="j")[['passed']])
    colors <- brewer.pal(3, "Dark2")[2:1]
    ggplot(dt, aes(value, fill=passed)) +
        geom_histogram(bins=bins) +
        scale_x_log10() +
        scale_y_log10() +
        scale_fill_manual(values=colors, name="Passed",
                labels=c("False", "True")) +
        xlab("Mean Junction Expression") +
        ylab("Count") +
        ggtitle("Expression filtering") +
        theme_bw() +
        theme(legend.position=legend.position)
}


#'
#' Plot count correlation
#'
#' Count correlation heatmap function
#'
#' @rdname plotFunctions
#' @export
plotCountCorHeatmap <- function(fds, type=c("psi5", "psi3", "psiSite"),
                    logit=FALSE, topN=50000, topJ=5000, minMedian=1,
                    main=NULL, normalized=FALSE, show_rownames=FALSE,
                    show_colnames=FALSE, minDeltaPsi=0.1, annotation_col=NA,
                    annotation_row=NA, border_color=NA, nClust=5,
                    plotType=c("sampleCorrelation", "junctionSample"),
                    sampleClustering=NULL, plotMeanPsi=TRUE, plotCov=TRUE, ...){

    type <- match.arg(type)
    plotType <- match.arg(plotType)

    # use counts as matrix, otherwise x(fds,...) does not work later on
    counts(fds, type=type, side="other", HDF5=FALSE)      <-
        as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <-
        as.matrix(counts(fds, type=type, side="ofInterest"))

    kmat <- as.matrix(counts(fds, type=type, side="ofIn"))
    nmat <- kmat + as.matrix(counts(fds, type=type, side="other"))

    expRowsMedian <- rowMedians(kmat) >= minMedian
    expRowsMax    <- rowMax(kmat)     >= 10
    table(expRowsMax & expRowsMedian)

    skmat <- kmat[expRowsMax & expRowsMedian,]
    snmat <- nmat[expRowsMax & expRowsMedian,]

    xmat <- (skmat + 1)/(snmat + 2)
    if(isTRUE(logit)){
        xmat <- qlogisWithCap(xmat)
    }
    xmat_rc    <- xmat - rowMeans(xmat)

    xmat_rc_sd <- rowSds(xmat_rc)
    plotIdx <- rank(xmat_rc_sd) >= length(xmat_rc_sd) - topN
    xmat_rc_2_plot <- xmat_rc[plotIdx,]
    cormatS <- cor(xmat_rc_2_plot)
    if(isTRUE(normalized)){
        pred_mu <- as.matrix(predictedMeans(fds, type=type)[
            expRowsMax & expRowsMedian,][plotIdx,])
        pred_mu <- qlogisWithCap(pred_mu)
        lpred_mu_rc <- pred_mu - rowMeans(pred_mu)
        xmat_rc_2_plot <- xmat_rc_2_plot - lpred_mu_rc
    }
    cormat <- cor(xmat_rc_2_plot)
    breaks <- seq(-1, 1, length.out=50)

    if(plotType == "junctionSample"){

        if(isTRUE(normalized)){
            pred_mu <- as.matrix(predictedMeans(fds, type=type)[
                expRowsMax & expRowsMedian,])
            if(isTRUE(logit)){
                pred_mu <- qlogisWithCap(pred_mu)
            }
            lpred_mu_rc <- pred_mu - rowMeans(pred_mu)
            xmat_rc <- xmat_rc - lpred_mu_rc
        }

        fds <- fds[expRowsMax & expRowsMedian,,by=type]
        j2keepVa <- variableJunctions(fds, type, minDeltaPsi)
        j2keepDP <- rowQuantiles(kmat[expRowsMax & expRowsMedian,],
                                 probs=0.75) >= 10
        j2keep <- j2keepDP & j2keepVa
        xmat_rc_2_plot <- xmat_rc[j2keep,]
        mostVarKeep <- subsetKMostVariableJunctions(fds[j2keep,,by=type],
                                                    type, topJ)
        xmat_rc_2_plot <- xmat_rc_2_plot[mostVarKeep,]
        rownames(xmat_rc_2_plot) <- seq_len(nrow(xmat_rc_2_plot))
        breaks <- seq(-5, 5, length.out=50)

    }


    if(is.character(annotation_col)){
        annotation_col <- getColDataAsDFFactors(fds, annotation_col)
    }
    if(is.character(annotation_row)){
        annotation_row <- getColDataAsDFFactors(fds, annotation_row)
    }

    # annotate with sample clusters
    if(is.null(sampleClustering)){
        # annotate samples with clusters from sample correlation heatmap
        clusters <- as.factor(cutree(hclust(dist(cormatS)), k=nClust))
    } else if(!is.na(sampleClustering)){
        clusters <- sampleClustering
    }

    if(!isTRUE(is.na(sampleClustering))){
        if(!is.null(nrow(annotation_col))){
            annotation_col$sampleCluster <- clusters
        } else {
            annotation_col <- data.frame(sampleCluster=clusters)
        }
    }


    if(plotType == "junctionSample"){

        # annotate junctions with meanPsi and meanCoverage
        xmat <- xmat[j2keep,]
        xmat <- xmat[mostVarKeep,]
        meanPsi <- if(isTRUE(logit)){
            rowMeans(plogis(xmat))
        } else{
            rowMeans(xmat)
        }
        meanPsiBins <- cut(meanPsi, breaks = c(0, 0.33, 0.66, 1),
                           include.lowest=TRUE)
        if(isTRUE(plotMeanPsi)){
            if(!is.null(nrow(annotation_row))){
                annotation_row$meanPsi <- meanPsiBins
            } else{
                annotation_row <- data.frame(meanPsi=meanPsiBins)
            }
        }

        snmat <- snmat[j2keep,]
        snmat <- snmat[mostVarKeep,]
        meanCoverage <- rowMeans(snmat)
        cutpoints <- sort(unique(round(log10(meanCoverage))))
        if(max(cutpoints) < ceiling(log10(max(meanCoverage)))){
            cutpoints <- c(cutpoints, ceiling(log10(max(meanCoverage))))
        }
        meanCoverage <- cut(meanCoverage, breaks=10^(cutpoints),
                            include.lowest=TRUE)

        if(isTRUE(plotCov)){
            annotation_row$meanCoverage <- meanCoverage
        }
        if(isTRUE(nrow(annotation_row) > 0)){
            rownames(annotation_row) <- rownames(xmat_rc_2_plot)
        }
        cormat <- xmat_rc_2_plot
    }

    if(is.null(main)){
        main <- ifelse(normalized, "Normalized row-centered ", 
                       "Raw row-centered ")
        if(plotType == "sampleCorrelation"){
            if(isTRUE(logit)){
                main <- paste0(main, "Logit(PSI) correlation (", type, ")")
            } else {
                main <- paste0(main, "PSI correlation (", type, ")")
            }
        } else {
            if(isTRUE(logit)){
                main <- paste0(main, "Logit(PSI) data (", type, ", top ", topJ, 
                               ")")
            } else {
                main <- paste0(main, "PSI data (", type, ", top ", topJ, ")")
            }
        }
    }

    pheatmap(cormat, show_rownames=show_rownames, show_colnames=show_colnames,
             main=main, annotation_col=annotation_col, breaks=breaks,
             annotation_row=annotation_row, ..., border_color=border_color,
             color=colorRampPalette(colors=rev(brewer.pal(11, "RdBu")))(50)
    )
}

#'
#' helper function to get the annotation as data frame from the col data object
#'
#' @noRd
getColDataAsDFFactors <- function(fds, names){
    tmpDF <- data.frame(colData(fds)[,names])
    colnames(tmpDF) <- names
    for(i in names){
        if(any(is.na(tmpDF[, i]))){
            tmpDF[,i] <- as.factor(paste0("", tmpDF[,i]))
        }
        if(is.integer(tmpDF[,i]) && length(levels(as.factor(tmpDF[,i]))) <= 10){
            tmpDF[,i] <- as.factor(paste0("", tmpDF[,i]))
        }
    }
    rownames(tmpDF) <- rownames(colData(fds))
    return(tmpDF)
}

#' @noRd
qlogisWithCap <- function(x){
    ans <- qlogis(x)
    ans[is.infinite(ans)] <- NA
    rowm <- rowMaxs(ans, na.rm=TRUE)
    idx <- which(is.na(ans), arr.ind=TRUE)
    ans[idx] <- rowm[idx[,"row"]]
    return(ans)
}

#'
#' Helper to get nice Splice metric labels in ggplot
#'
#' @noRd
ggplotLabelPsi <- function(type){
    vapply(type, FUN=function(x)
        switch (x,
                psi5 = c(bquote(psi[5])),
                psi3 = c(bquote(psi[3])),
                psiSite = c(bquote(theta))),
        FUN.VALUE=c(bquote(psi[3])))
}

#' @noRd
plotLoss <- function(fds, type){
    type
    fds <- fds_new
    lossList <- metadata(fds)[[paste0('loss_', type)]]
    dt2plot <- as.data.table(melt(do.call(rbind, list(
        max=colMaxs(lossList, na.rm=TRUE),
        mean=colMeans(lossList, na.rm=TRUE),
        min=colMins(lossList, na.rm=TRUE),
        sd=colSds(lossList, na.rm=TRUE))), value.name="Loss"))
    colnames(dt2plot)[c(1,2)] <- c("Aggregation", "IterSteps")
    
    # set iterations
    dt2plot[grepl("init",IterSteps),iteration:=0]
    dt2plot[is.na(iteration),iteration:=as.numeric(as.character(
        gsub("_.*", "" ,gsub("final_", "" , IterSteps))))]
    
    # set step
    dt2plot[grepl("final", IterSteps),Step:="D fit"]
    dt2plot[is.na(Step) & !grepl("init", IterSteps),Step:="E fit"]
    dt2plot[is.na(Step),Step:="Init"]
    dt2plot[,Iteration:=seq_len(.N)/1,by="Aggregation,Step"]
    dt2plot[Step=="E fit", Iteration:=Iteration/3]
    dt2plot[Step=="D fit", Iteration:=Iteration/2]
    
    # summaries
    ggplot(dt2plot[Aggregation != "sd"], aes(x=Iteration, y=Loss, col=Step, 
                                             lty=Aggregation)) +
        geom_ribbon(data=dt2plot[Aggregation == "mean",], aes(
            ymin = dt2plot[Aggregation == "mean", Loss] - 
                dt2plot[Aggregation == "sd", Loss],
            ymax = dt2plot[Aggregation == "mean", Loss] + 
                dt2plot[Aggregation == "sd", Loss]), fill="gray80", col=NA) +
        geom_line() +
        scale_y_log10()
}
