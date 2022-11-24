#' 
#' Visualization functions for FRASER
#'
#' The FRASER package provides mutliple functions to visualize
#' the data and the results of a full data set analysis.
#'
#' This is the list of all plotting function provided by FRASER:
#' \itemize{
#'   \item plotAberrantPerSample()
#'   \item plotVolcano()
#'   \item plotExpression()
#'   \item plotQQ()
#'   \item plotExpectedVsObservedPsi()
#'   \item plotCountCorHeatmap()
#'   \item plotFilterExpression()
#'   \item plotFilterVariability()
#'   \item plotEncDimSearch()
#'   \item plotBamCoverage()
#'   \item plotBamCoverageFromResultTable()
#'   \item plotManhattan()
#' }
#'
#' For a detailed description of each plot function please see the details.
#' Most of the functions share the same parameters.
#'
#### Data specific parameters
#' @param object,fds An \code{\link{FraserDataSet}} object.
#' @param type The psi type: either psi5, psi3 or theta (for SE).
#' @param sampleID A sample ID which should be plotted. Can also be a vector.
#'             Integers are treated as indices.
#' @param idx A junction site ID or gene ID or one of both, which
#'             should be plotted. Can also be a vector. Integers are treated
#'             as indices.
#' @param padjCutoff,zScoreCutoff,deltaPsiCutoff Significance, Z-score or delta
#'             psi cutoff to mark outliers
#' @param global Flag to plot a global Q-Q plot, default FALSE
#' @param normalized If TRUE, the normalized psi values are used, the default,
#'             otherwise the raw psi values
#' @param aggregate If TRUE, the pvalues are aggregated by gene (default), 
#'             otherwise junction level pvalues are used (default for Q-Q plot).
#' @param result The result table to be used by the method.
#' @param label Indicates the genes or samples that will be labelled in the 
#'             plot (only for \code{basePlot=TRUE}). Setting 
#'             \code{label="aberrant"} will label all aberrant genes or 
#'             samples. Labelling can be turned off by setting 
#'             \code{label=NULL}. The user can also provide a custom 
#'             list of gene symbols or sampleIDs.
#' @param BPPARAM BiocParallel parameter to use.
#' @param Ncpus Number of cores to use.
#' @param plotType The type of plot that should be shown as character string. 
#'             For plotEncDimSearch, it has to be either \code{"auc"} 
#'             for a plot of the area under the curve (AUC) or 
#'             \code{"loss"} for the model loss. For the correlation heatmap,
#'             it can be either \code{"sampleCorrelation"} for a
#'             sample-sample correlation heatmap or \code{"junctionSample"}
#'             for a junction-sample correlation heatmap.
#' @param onlyVariableIntrons Logical value indicating whether to show only 
#'             introns that also pass the variability filter. Defaults to 
#'             FALSE.
#' @param onlyExpressedIntrons Logical value indicating whether to show only 
#'             introns that also pass the expression filter. Defaults to 
#'             FALSE.
#' @param gr A GRanges object indicating the genomic range that should be shown 
#'             in \code{plotBamCoverage}.
#' @param control_samples The sampleIDs of the samples used as control in  
#'             \code{plotBamCoverage}.
#' @param min_junction_count The minimal junction count across samples required 
#'             for a junction to appear in the splicegraph and coverage tracks 
#'             of \code{plotBamCoverage}.
#' @param txdb A TxDb object giving the gene/transcript annotation to use.
#' @param orgDb A OrgDb object giving the mapping of gene ids and symbols.
#' @param show_full_gene Should the full genomic range of the gene be shown in 
#'             \code{plotBamCoverageFromResultTable} (default: FALSE)? 
#'             If FALSE, only a certain region (see parameters left_extension 
#'             and right_extension) around the outlier junction is shown. 
#' @param left_extension Indicating how far the plotted range around the outlier 
#'             junction should be extended to the left in 
#'             \code{plotBamCoverageFromResultTable}.
#' @param right_extension Indicating how far the plotted range around the 
#'             outlier junction should be extended to the right in 
#'             \code{plotBamCoverageFromResultTable}.
#' @param res_gene_col The column name in the given results table that 
#'             contains the gene annotation.
#' @param res_geneid_type The type of gene annotation in the results table in 
#'             \code{res_gene_col} (e.g. SYMBOL or ENTREZID etc.). This 
#'             information is needed for mapping between the results table and 
#'             the provided annotation in the txdb object. 
#' @param txdb_geneid_type The type of gene_id present in \code{genes(txdb)} 
#'             (e.g. ENTREZID). This information is needed for 
#'             mapping between the results table and the provided annotation 
#'             in the txdb object. 
#' @param highlight_range A \code{GenomicRanges} or \code{GenomicRangesList} 
#'             object of ranges to be highlighted in the splicegraph of 
#'             \code{plotBamCoverage}.
#' @param highlight_range_color The color of highlighted ranges in 
#'             the splicegraph of \code{plotBamCoverage}.
#' @param toscale In \code{plotBamCoverage}, indicates which part of the 
#'             plotted region should be drawn to scale. Possible values are 
#'             'exon' (exonic regions are drawn to scale), 
#'             'gene' (both exonic and intronic regions are drawn to scale) or 
#'             'none' (exonic and intronic regions have constant length) 
#'             (see SGSeq package).
#' @param splicegraph_labels Indicated the format of exon/splice junction 
#'             labels in the splicegraph of \code{plotBamCoverage}. 
#'             Possible values are 'genomic_range' (gives the start position 
#'             of the first exon and the end position of the last exon that 
#'             are shown),  'id' (format E1,... J1,...), 'name' (format 
#'             type:chromosome:start-end:strand for each feature), 
#'             'none' for no labels (see SGSeq package).             
#' @param splicegraph_position The position of the splicegraph relative to the 
#'             coverage tracks in \code{plotBamCoverage}. Possible values 
#'             are 'top' (default) and 'bottom'.
#'              
#### Graphical parameters
#' @param main Title for the plot, if missing a default title will be used.
#' @param colGroup Group of samples that should be colored.
#' @param basePlot if \code{TRUE} (default), use the R base plot version, else
#'             use the plotly framework.
#' @param conf.alpha If set, a confidence interval is plotted, defaults to 0.05
#' @param samplingPrecision Plot only non overlapping points in Q-Q plot to 
#'             reduce number of points to plot. Defines the digits to round to. 
#' @param logit If TRUE, the default, psi values are plotted in logit space.
#' @param nClust Number of clusters to show in the row and
#'             column dendrograms.
#' @param sampleClustering A clustering of the samples that should be used as an
#'             annotation of the heatmap.
#' @param show_rownames,show_colnames Logical value indicating whether to show
#'             row or column names on the heatmap axes.
#' @param annotation_col,annotation_row Row or column annotations that should be
#'             plotted on the heatmap.
#' @param topN Top x most variable junctions that should be used for the
#'             calculation of sample x sample correlations. 
#' @param topJ Top x most variable junctions that should be displayed in the
#'             junction-sample correlation heatmap. Only applies if plotType 
#'             is "junctionSample".
#' @param minMedian,minCount,minDeltaPsi Minimal median (\eqn{m \ge 1}), 
#'             delta psi (\eqn{|\Delta\psi| > 0.1}), read count (\eqn{n \ge 10})
#'             value of a junction to be considered for the correlation heatmap.
#' @param border_color Sets the border color of the heatmap
#' @param plotMeanPsi,plotCov If \code{TRUE}, then the heatmap is annotated with
#'             the mean psi values or the junction coverage.
#' @param bins Set the number of bins to be used in the histogram.
#' @param legend.position Set legend position (x and y coordinate), defaults to
#'             the top right corner.
#' @param color_annotated The color for exons and junctions present in 
#'             the given annotation (in the splicegraph of 
#'             \code{plotBamCoverage}).
#' @param color_novel The color for novel exons and junctions not present in 
#'             the given annotation (in the splicegraph of 
#'             \code{plotBamCoverage}).
#' @param color_sample_interest The color in \code{plotBamCoverage} for the 
#'             sample of interest.
#' @param color_control_samples The color in \code{plotBamCoverage} for the 
#'             samples used as controls.
#' @param curvature_splicegraph The curvature of the junction arcs in the 
#'             splicegraph in \code{plotBamCoverage}. Decrease this value 
#'             for flatter arcs and increase it for steeper arcs.
#' @param curvature_coverage The curvature of the junction arcs in the 
#'             coverage tracks of \code{plotBamCoverage}. Decrease this 
#'             value for flatter arcs and increase it for steeper arcs.
#' @param mar The margin of the plot area for \code{plotBamCoverage} 
#'             (b,l,t,r).
#' @param cex For controlling the size of text and numbers in 
#'             \code{plotBamCoverage}.
#' @param color_chr Interchanging colors by chromosome for \code{plotManhattan}.
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
#' given psi type (\eqn{\psi_5, \psi_3, or \theta}{\psi5, \psi3, or \theta}
#' (SE)) for all samples.
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
#' FraserDataSet object contains the \code{passedFilter} column, it will plot 
#' both FPKM distributions for the expressed introns and for the filtered 
#' introns.
#' 
#' \code{plotFilterVariability}: The distribution of maximal delta Psi values. 
#' If the FraserDataSet object contains the \code{passedFilter} column, 
#' it will plot both maximal delta Psi distributions for the variable 
#' introns and for the filtered (i.e. non-variable) introns.
#'
#' \code{plotEncDimSearch}: Visualization of the hyperparameter optimization.
#' It plots the encoding dimension against the achieved loss (area under the
#' precision-recall curve). From this plot the optimum should be choosen for
#' the \code{q} in fitting process.
#' 
#' \code{plotManhattan}: A Manhattan plot showing the junction pvalues by 
#' genomic position. Useful to identify if outliers cluster by genomic position.
#' 
#' \code{plotBamCoverage}: A sashimi plot showing the read coverage from 
#' the underlying bam files for a given genomic range and sampleIDs. 
#' 
#' \code{plotBamCoverageFromResultTable}: A sashimi plot showing the read 
#' coverage from the underlying bam files for a row in the results table. Can 
#' either show the full range of the gene with the outlier junction or only a 
#' certain region around the outlier. 
#'
#' @return If base R graphics are used nothing is returned else the plotly or
#'             the gplot object is returned.
#'
#' @name plotFunctions
#' @rdname plotFunctions
#' @aliases plotFunctions plotAberrantPerSample plotVolcano plotQQ 
#'             plotExpression plotCountCorHeatmap plotFilterExpression 
#'             plotExpectedVsObservedPsi plotEncDimSearch plotManhattan
#'             plotBamCoverage plotBamCoverageFromResultTable
#' @examples
#' # create full FRASER object 
#' fds <- makeSimulatedFraserDataSet(m=40, j=200)
#' fds <- calculatePSIValues(fds)
#' fds <- filterExpressionAndVariability(fds, filter=FALSE)
#' # this step should be done for more dimensions in practice
#' fds <- optimHyperParams(fds, "jaccard", q_param=c(2,5,10,25))
#' fds <- FRASER(fds)
#' 
#' # QC plotting
#' plotFilterExpression(fds)
#' plotFilterVariability(fds)
#' plotCountCorHeatmap(fds, "jaccard")
#' plotCountCorHeatmap(fds, "jaccard", normalized=TRUE)
#' plotEncDimSearch(fds, type="jaccard")
#' 
#' # extract results 
#' plotAberrantPerSample(fds, aggregate=FALSE)
#' plotVolcano(fds, "sample1", "jaccard")
#' 
#' # dive into gene/sample level results
#' res <- results(fds)
#' res
#' plotExpression(fds, result=res[1])
#' plotQQ(fds, result=res[1])
#' plotExpectedVsObservedPsi(fds, res=res[1])
#' 
#' # create manhattan plot of pvalues by genomic position
#' plotManhattan(fds, type="jaccard", sampleID="sample10")
#' 
#' # plot splice graph and coverage from bam files in a given region
#' if(require(SGSeq)){
#'     fds <- createTestFraserSettings()
#'     gr <- GRanges(seqnames="chr19", 
#'         IRanges(start=7587496, end=7598895), 
#'         strand="+")
#'     plotBamCoverage(fds, gr=gr, sampleID="sample3", 
#'         control_samples="sample2", min_junction_count=5,
#'         curvature_splicegraph=1, curvature_coverage=1, 
#'         mar=c(1, 7, 0.1, 3))
#' 
#'     # plot coverage from bam file for a row in the result table
#'     fds <- createTestFraserDataSet()
#'     require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'     require(org.Hs.eg.db)
#'     orgDb <- org.Hs.eg.db
#'  
#'     res <- results(fds, padjCutoff=NA, deltaPsiCutoff=NA, zScoreCutoff=NA)
#'     res_dt <- as.data.table(res)
#'     res_dt <- res_dt[sampleID == "sample2",]
#'     
#'     # plot full range of gene containing outlier junction
#'     plotBamCoverageFromResultTable(fds, result=res_dt[1,], show_full_gene=TRUE,
#'         txdb=txdb, orgDb=orgDb, control_samples="sample3")
#'     
#'     # plot only certain range around outlier junction
#'     plotBamCoverageFromResultTable(fds, result=res_dt[1,], show_full_gene=FALSE, 
#'         control_samples="sample3", curvature_splicegraph=0.5, txdb=txdb,
#'         curvature_coverage=0.5, right_extension=5000, left_extension=5000,
#'         splicegraph_labels="id")
#' }
#' 
NULL


plotVolcano.FRASER <- function(object, sampleID, 
                    type=psiTypes_avail, basePlot=TRUE, 
                    aggregate=FALSE, main=NULL, label=NULL,
                    deltaPsiCutoff=0.1, padjCutoff=0.1, ...){
    
    type <- match.arg(type)

    dt <- getPlottingDT(object, axis="col", type=type, idx=sampleID,
            aggregate=aggregate, deltaPsiCutoff=deltaPsiCutoff, 
            padjCutoff=padjCutoff, ...)
    
    g <- ggplot(dt, aes(x=deltaPsi, y=-log10(pval), color=aberrant, 
                        label=featureID, text=paste0(
            "SampleID: ", sampleID, "<br>",
            "featureID: ", featureID, "<br>",
            "Raw count (K): ", k, "<br>",
            "Raw total count (N): ", n, "<br>",
            "p value: ", signif(pval, 5), "<br>",
            "delta Psi: ", round(deltaPsi, 2), "<br>",
            "Type: ", type))) +
        geom_point(aes(alpha=ifelse(aberrant == TRUE, 1, 0.8))) +
        xlab(as.expression(
                bquote(paste(Delta, .(ggplotLabelPsi(type)[[1]]) ))
            )) +
        ylab(expression(paste(-log[10], "(P value)"))) +
        theme_bw() +
        theme(legend.position="none") +
        scale_color_manual(values=c("gray40", "firebrick"))
    
    if(!is.na(deltaPsiCutoff)){
        g <- g + 
            geom_vline(xintercept=c(-deltaPsiCutoff, deltaPsiCutoff),
                    color="firebrick", linetype=2)
    }
    
    if(!is.na(padjCutoff)){
        if(dt[padj <= padjCutoff, .N] > 0){
            padj_line <- min(dt[padj <= padjCutoff, -log10(pval)])
            padj_line <- min(dt[padj <= padjCutoff, -log10(pval)])
        }
        if(!"padj_line" %in% ls() || padj_line > 10 || is.na(padj_line)){
            padj_line <- 6
        }
        g <- g + 
            geom_hline(yintercept=padj_line, color="firebrick", linetype=4)
    }
    
    
    if(isFALSE(basePlot)){
        g <- g + xlab(paste("delta", 
                            ggplotLabelPsi(type, asCharacter=TRUE)[[1]])) +
            ylab("-log[10](P value)")
        if(is.null(main)){
            main <- paste0("Volcano plot: ", sampleID, ", ", 
                            ggplotLabelPsi(type, asCharacter=TRUE)[[1]])
        }
    } else{
        if(!is.null(label)){
            if(isScalarCharacter(label) && label == "aberrant"){
                if(nrow(dt[aberrant == TRUE,]) > 0){
                    g <- g + geom_text_repel(data=dt[aberrant == TRUE,],
                                        fontface='bold', hjust=-.2, vjust=-.2)
                }
            }
            else{
                if(nrow(dt[featureID %in% label]) > 0){
                    g <- g + geom_text_repel(data=
                                        subset(dt, featureID %in% label),
                                    fontface='bold', hjust=-.2, vjust=-.2)
                }
                if(any(!(label %in% dt[,featureID]))){
                    warning("Did not find gene(s) ", 
                            paste(label[!(label %in% dt[,featureID])], 
                                collapse=", "), " to label.")
                }
            }
        }
        
        if(is.null(main)){
            main <- as.expression(bquote(paste(
                bold("Volcano plot: "), .(sampleID), ", ",
                .(ggplotLabelPsi(type)[[1]]))))
        }  
    }
    g <- g + ggtitle(main)
    
    plotBasePlot(g, basePlot)
}

#'
#' Volcano plot
#'
#' Plots the p values over the delta psi values, known as volcano plot.
#' Visualizes per sample the outliers. By type and aggregate by
#' gene if requested.
#'
#' @rdname plotFunctions
#' @export
setMethod("plotVolcano", signature="FraserDataSet", plotVolcano.FRASER)


plotAberrantPerSample.FRASER <- function(object, main, 
                    type=psiTypes_avail,
                    padjCutoff=0.1, zScoreCutoff=NA, deltaPsiCutoff=0.1,
                    aggregate=TRUE, BPPARAM=bpparam(), ...){

    type <- match.arg(type, several.ok=TRUE)

    if(missing(main)){
        main <- paste('Aberrant events per sample')
        if(isTRUE(aggregate)){
            main <- paste(main, "by gene")
        }
    }

    # extract outliers
    outliers <- bplapply(type, aberrant, object=object, by="sample",
            padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
            deltaPsiCutoff=deltaPsiCutoff, aggregate=aggregate, ..., 
            BPPARAM=BPPARAM)
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
        scale_color_brewer(palette="Dark2", name=element_blank(),
                labels=ggplotLabelPsi(dt2p[,unique(type)])) +
        scale_linetype_manual(name="", values=2, labels="Median")

    if(!all(dt2p[,value] == 0)){
        g <- g + scale_y_log10()
    }
    
    g
}

#'
#' Number of aberrant events per sample
#'
#' Plot the number of aberrant events per samples
#'
#' @rdname plotFunctions
#' @export
setMethod("plotAberrantPerSample", signature="FraserDataSet",
        plotAberrantPerSample.FRASER)


#'
#' Junction expression plot
#'
#' Plots the observed split reads of the junction of interest over all reads
#' coming from the given donor/acceptor.
#'
#' @rdname plotFunctions
#' @export
plotExpression <- function(fds, type=psiTypes_avail,
                    idx=NULL, result=NULL, colGroup=NULL, 
                    basePlot=TRUE, main=NULL, label="aberrant", ...){
    if(!is.null(result)){
        type <- as.character(result$type)
        idx <- getIndexFromResultTable(fds, result)
    } else {
        type <- match.arg(type)
    }

    dt <- getPlottingDT(fds, axis="row", type=type, idx=idx, ...)
    dt[,featureID:=limitGeneNamesList(featureID, maxLength=3)]
    
    if(!is.null(colGroup)){
        if(all(colGroup %in% samples(fds))){
            colGroup <- samples(fds) %in% colGroup
        }
        dt[colGroup,aberrant:=TRUE]
    }
    dt[,aberrant:=factor(aberrant, levels=c("TRUE", "FALSE"))]

    if(is.null(main)){
        if(isTRUE(basePlot)){
            main <- as.expression(bquote(bold(paste(
                .(ggplotLabelPsi(type)[[1]]), " expression plot: ",
                bolditalic(.(as.character(dt[,unique(featureID)]))),
                " (site ", .(as.character(dt[,unique(idx)])), ")"))))
        } else{
            main <- paste0(ggplotLabelPsi(type, asCharacter=TRUE)[[1]], 
                        " expression plot: ", dt[,unique(featureID)], 
                        " (site ", dt[,unique(idx)], ")")
        }
    }

    g <- ggplot(dt, aes(x=n + 2, y=k + 1, color=aberrant, label=sampleID, 
                        text=paste0(
            "Sample: ", sampleID, "<br>",
            "Counts (K): ", k, "<br>",
            "Total counts (N): ", n, "<br>",
            "p value: ", signif(pval, 5), "<br>",
            "padjust: ", signif(padj, 5), "<br>",
            "Observed Psi: ", round(obsPsi, 2), "<br>",
            "Predicted mu: ", round(predPsi, 2), "<br>"))) +
        geom_point(alpha=ifelse(as.character(dt$aberrant) == "TRUE", 1, 0.7)) +
        scale_x_log10() +
        scale_y_log10() +
        theme_bw() +
        theme(legend.position="none", title=) +
        xlab("Total junction coverage + 2 (N)") +
        ylab("Junction count + 1 (K)") +
        ggtitle(main) +
        annotation_logticks(sides='bl')
    
    if(isTRUE(basePlot) && !is.null(label)){
        if(isScalarCharacter(label) && label == "aberrant"){
            if(nrow(dt[aberrant == TRUE,]) > 0){
                g <- g + geom_text_repel(data=dt[aberrant == TRUE,], 
                                        aes(col=aberrant),
                                        fontface='bold', hjust=-.2, vjust=-.2)
            }
        }
        else{
            if(nrow(dt[sampleID %in% label]) > 0){
                g <- g + geom_text_repel(data=subset(dt, sampleID %in% label), 
                                        aes(col=aberrant),
                                        fontface='bold', hjust=-.2, vjust=-.2)
            }
            if(any(!(label %in% dt[,sampleID]))){
                warning("Did not find sample(s) ", 
                        paste(label[!(label %in% dt[,sampleID])], 
                            collapse=", "), " to label.")
            }
        }
    }

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
plotExpectedVsObservedPsi <- function(fds, type=psiTypes_avail,
                    idx=NULL, result=NULL, colGroup=NULL, main=NULL,
                    basePlot=TRUE, label="aberrant", ...){
    type <- match.arg(type)
    
    # get plotting data
    dt   <- getPlottingDT(fds, axis="row", type=type, result=result, 
            idx=idx, ...)
    type <- as.character(unique(dt$type))
    idx  <- unique(dt$idx)
    dt[,featureID:=limitGeneNamesList(featureID, maxLength=3)]
    
    if(is.null(main)){
        if(isTRUE(basePlot)){
            main <- as.expression(bquote(bold(paste(
                .(ggplotLabelPsi(type)[[1]]), 
                " observed expression vs prediction plot: ",
                bolditalic(.(as.character(dt[,unique(featureID)]))),
                " (site ", .(as.character(idx)), ")"))))
        } else{
            main <- paste0(ggplotLabelPsi(type, asCharacter=TRUE)[[1]], 
                            " observed expression vs prediction plot: ", 
                            dt[,unique(featureID)], " (site ", idx, ")")
        }
    }
    
    if(!is.null(colGroup)){
        if(is.logical(colGroup)){
            dt[colGroup, aberrant:=TRUE]
        } else {
            warning("not implemented yet!")
        }
    }

    if(isTRUE(basePlot)){
        ylab <- bquote("Observed " ~ .(ggplotLabelPsi(type)[[1]]))
        xlab <- bquote("Predicted " ~ .(ggplotLabelPsi(type)[[1]]))
    } else{
        ylab <- paste("Observed", ggplotLabelPsi(type, asCharacter=TRUE)[[1]])
        xlab <- paste("Predicted", ggplotLabelPsi(type, asCharacter=TRUE)[[1]])
    }
    
    g <- ggplot(dt, aes(y=obsPsi, x=predPsi, label=sampleID, color=aberrant, 
            text=paste0(
            "Sample: ", sampleID, "<br>",
            "Counts (K): ", k, "<br>",
            "Total counts (N): ", n, "<br>",
            "p value: ", signif(pval, 5), "<br>",
            "padjust: ", signif(padj, 5), "<br>",
            "Observed Psi: ", round(obsPsi, 2), "<br>",
            "Predicted mu: ", round(predPsi, 2), "<br>"))) +
        geom_point(alpha=ifelse(dt$aberrant, 1, 0.5),
                color=c("gray70", "firebrick")[dt$aberrant + 1]) +
        geom_abline(intercept = 0, slope=1) +
        xlim(c(0,1)) + ylim(c(0,1)) +
        theme_bw() +
        theme(legend.position="none") +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(main)
    
    if(isTRUE(basePlot) && !is.null(label)){
        if(isScalarCharacter(label) && label == "aberrant"){
            if(nrow(dt[aberrant == TRUE,]) > 0){
                g <- g + geom_text_repel(data=dt[aberrant == TRUE,], 
                                        aes(col=aberrant),
                                        fontface='bold', hjust=-.2, vjust=-.2)
            }
        }
        else{
            if(nrow(dt[sampleID %in% label]) > 0){
                g <- g + geom_text_repel(data=subset(dt, sampleID %in% label), 
                                    aes(col=aberrant),
                                    fontface='bold', hjust=-.2, vjust=-.2)
            }
            if(any(!(label %in% dt[,sampleID]))){
                warning("Did not find sample(s) ", 
                        paste(label[!(label %in% dt[,sampleID])], 
                                collapse=", "), " to label.")
            }
        }
    }    

    if(is.null(colGroup)){
        g <- g + scale_colour_manual(
            values=c("FALSE"="gray70", "TRUE"="firebrick"))
    }
    
    plotBasePlot(g, basePlot)
}


plotQQ.FRASER <- function(object, type=NULL, idx=NULL, result=NULL, 
                    aggregate=FALSE, global=FALSE, main=NULL, conf.alpha=0.05,
                    samplingPrecision=3, basePlot=TRUE, label="aberrant",
                    Ncpus=min(3, getDTthreads()), ...){

    # check parameters
    if(is.null(aggregate)){
        aggregate <- isTRUE(global)
    } else if(!(is.logical(aggregate) |
                all(aggregate %in% colnames(mcols(object))))){
        stop("Please provide TRUE/FALSE or a ",
            "charactor matching a column name in mcols.")
    }

    if(isTRUE(global)){
        if(is.null(type)){
            type <- psiTypes
        }
        dt <- rbindlist(bplapply(type, getPlottingDT, fds=object, axis="col",
                idx=TRUE, aggregate=aggregate, Ncpus=Ncpus, ...))
        # remove duplicated entries donor/acceptor sites if not aggregated 
        # by a feature
        if(isFALSE(aggregate)){
            dt <- dt[!duplicated(dt, by=c("type", "spliceID", "sampleID"))]
        }
    } else {
        dots <- list(...)
        if(!"pvalLevel" %in% names(dots)){
            dots[["pvalLevel"]] <- "junction"
        }
        dots <- append(list(fds=object, axis="row", type=type, idx=idx, 
                            result=result, aggregate=aggregate), 
                        dots)
        dt <- do.call(getPlottingDT, args=dots)
    }
    if(is.null(main)){
        if(isTRUE(global)){
            main <- "Global QQ plot"
        } else {
            type <- as.character(dt[,unique(type)])
            featureID <- as.character(dt[,unique(featureID)])
            featureID <- limitGeneNamesList(featureID, maxLength=3)
            if(isTRUE(basePlot)){
                main <- as.expression(bquote(bold(paste(
                        .(ggplotLabelPsi(type)[[1]]),
                        " Q-Q plot: ", bolditalic(.(featureID)),
                        " (site ", .(as.character(dt[,unique(idx)])), ")"))))
            } else{
                main <- paste0(ggplotLabelPsi(type, asCharacter=TRUE)[[1]],
                                " Q-Q plot: ", featureID, 
                                " (site ", dt[,unique(idx)], ")")
            }
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
    if(isTRUE(samplingPrecision) | isScalarNumeric(samplingPrecision)){
        if(is.logical(samplingPrecision)){
            samplingPrecision <- 3
        }
        mypoints <- !duplicated(dt2p[,.(obs=round(obs, samplingPrecision), 
                        exp=round(exp, samplingPrecision), type)], 
                by=c("obs", "exp", "type"))
        dt2p[,plotPoint:=mypoints]
    }

    # create qq-plot
    g <- ggplot(dt2p[plotPoint == TRUE,], aes(x=exp, y=obs, col=aberrant, 
            label=sampleID,
            text=paste(
                "<br>SampleID: ", sampleID, "<br>K: ", k, "<br>N: ", n))) +
        geom_point() +
        theme_bw() +
        ggtitle(main) 
    
    if(isTRUE(basePlot)){
        g <- g +
            xlab(expression(-log[10]~"(expected P)")) +
            ylab(expression(-log[10]~"(observed P)"))
    } else{
        g <- g +
            xlab("-log[10] (expected P)") +
            ylab("-log[10] (observed P)")
    }

    # Set color scale for global/local
    if(isFALSE(global)){
        g <- g + scale_color_manual(values=c("black", "firebrick"),
                name="Aberrant") +
            theme(legend.position="none")
    } else {
        g$mapping$colour <- quote(type)
        g <- g + scale_color_brewer(palette="Dark2", name="Splice metric",
                labels=ggplotLabelPsi(unique(dt2p$type)))
    }


    # add confidence band if requesded
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(is.numeric(conf.alpha)){
        dt2p[,rank:=seq_len(.N), by=type]
        dt2p[plotPoint == TRUE,upper:=-log10(
                qbeta(conf.alpha/2, rank, max(rank) - rank)), by=type]
        dt2p[plotPoint == TRUE,lower:=-log10(
                qbeta(1-conf.alpha/2, rank, max(rank) - rank)), by=type]
        # only plot one psiType if multiple are existing
        if(length(unique(dt2p$type)) > 1){
            typeOrder <- c("theta", "psi5", "psi3")
            type2take <- min(which(typeOrder %in% unique(dt2p$type)))
            dt2p[type != typeOrder[type2take], 
                    c("upper", "lower"):=list(NA, NA)]
        }
        g <- g + geom_ribbon(data=dt2p[plotPoint == TRUE & !is.na(upper),],
                aes(x=exp, ymin=lower, ymax=upper, text=NULL),
                alpha=0.2, color="gray")
    }
    
    # add labels if requested
    if(isFALSE(global) && isTRUE(basePlot) && !is.null(label)){
        if(isScalarCharacter(label) && label == "aberrant"){
            if(nrow(dt2p[aberrant == TRUE,]) > 0){
                g <- g + geom_text_repel(data=dt2p[aberrant == TRUE,], 
                                        aes(col=aberrant),
                                        fontface='bold', hjust=-.2, vjust=-.2)
            }
        }
        else{
            if(nrow(dt2p[sampleID %in% label]) > 0){
                g <- g + geom_text_repel(data=subset(dt2p, sampleID %in% 
                                                                        label), 
                                        aes(col=aberrant),
                                        fontface='bold', hjust=-.2, vjust=-.2)
            }
            if(any(!(label %in% dt2p[,sampleID]))){
                warning("Did not find sample(s) ", 
                        paste(label[!(label %in% dt2p[,sampleID])], 
                                collapse=", "), " to label.")
            }
        }
    }   

    # add abline in the end
    g <- g + geom_abline(intercept=0, slope=1, col="firebrick")
    if(isFALSE(global)){
        return(plotBasePlot(g, basePlot))
    }
    g
}

#'
#' Q-Q plot
#'
#' Plots the quantile-quantile plot
#'
#' @rdname plotFunctions
#' @export
setMethod("plotQQ", signature="FraserDataSet", plotQQ.FRASER)


plotEncDimSearch.FRASER <- function(object, type=psiTypes_avail, 
                    plotType=c("auc", "loss")){
    type <- match.arg(type)
    plotType <- match.arg(plotType)
    data <- hyperParams(object, type=type, all=TRUE)
    if (is.null(data)) {
        warning(paste("no hyperparameters were estimated for", type, 
                        "\nPlease use `optimHyperParams` to compute them."))
        return(NULL)
    }
    if(!"nsubset" %in% colnames(data)){
        data[,nsubset:="NA"]
    }
    data[,noise:=as.factor(noise)]
    data[,nsubset:=as.factor(nsubset)]
    
    data[,isOptimalQ:= q == .SD[aroc == max(aroc), q], by="nsubset,noise"]

    if(plotType == "auc"){
        
        g1 <- ggplot(data, aes(q, aroc, col=nsubset, linetype=noise)) +
            geom_point() +
            geom_smooth(method="loess", formula=y~x) +
            geom_vline(data=data[isOptimalQ == TRUE,], 
                    mapping=aes(xintercept=q, col=nsubset, linetype=noise)) +
            ggtitle(as.expression(bquote(bold(paste(
                "Q estimation for ", .(ggplotLabelPsi(type)[[1]])))))) +
            xlab("Estimated q") +
            ylab("Area under the PR curve") +
            theme_bw(base_size=16)
        
        if(data[,uniqueN(nsubset)] == 1 & data[,uniqueN(noise)] == 1){
            g1 <- g1 + theme(legend.position='none')
        }
        
        g1
    }
    else{

        g2 <- ggplot(data, aes(q, eval, col=nsubset, linetype=noise)) +
            geom_point() +
            geom_smooth() +
            geom_vline(data=data[isOptimalQ == TRUE,], 
                    mapping=aes(xintercept=q, col=nsubset, linetype=noise)) +
            ggtitle(as.expression(bquote(bold(paste(
                "Q estimation for ", .(ggplotLabelPsi(type)[[1]])))))) +
            xlab("Estimated q") +
            ylab("Model loss") +
            theme_bw(base_size=16)
        
        if(data[,uniqueN(nsubset)] == 1 & data[,uniqueN(noise)] == 1){
            g2 <- g2 + theme(legend.position='none')
        }
        
        g2
    }

}

#'
#' Plots the results from the hyperparamter optimization.
#'
#' @rdname plotFunctions
#' @export
setMethod("plotEncDimSearch", signature="FraserDataSet", 
        plotEncDimSearch.FRASER)


#'
#' Plot filter expression
#'
#' Histogram of the geometric mean per junction based on the filter status
#'
#' @rdname plotFunctions
#' @export
plotFilterExpression <- function(fds, bins=200, legend.position=c(0.8, 0.8),
                                    onlyVariableIntrons=FALSE){
    
    # check that expression filter has been calculated
    if(!("passedExpression" %in% colnames(mcols(fds, type="j")))){
        stop("Please calculate the expression filter values first with the ",
                "filterExpression function.")
    }
    
    # get mean count for all junctions in the fds object
    cts    <- K(fds, "psi5")
    rowlgm <- exp(rowMeans(log(cts + 1)))

    dt <- data.table(
            value=rowlgm,
            passed=mcols(fds, type="j")[['passedExpression']])
    if(isTRUE(onlyVariableIntrons)){
        dt[,passed:=mcols(fds, type="j")[['passed']]]
    }
    dt[,passed:=factor(passed, levels=c(TRUE, FALSE))]
    
    colors <- brewer.pal(3, "Dark2")[seq_len(2)]
    ggplot(dt, aes(value, fill=passed)) +
        geom_histogram(bins=bins) +
        scale_x_log10() +
        scale_y_log10() +
        scale_fill_manual(values=colors, name="Passed",
                labels=c("True", "False")) +
        xlab("Mean Junction Expression") +
        ylab("Count") +
        ggtitle("Expression filtering") +
        theme_bw() +
        theme(legend.position=legend.position)
}

#'
#' Plot filter variability
#'
#' Histogram of minimal delta psi per junction
#'
#' @rdname plotFunctions
#' @export
plotFilterVariability <- function(fds, bins=200, legend.position=c(0.8, 0.8),
                                    onlyExpressedIntrons=FALSE){
    
    # check that expression filter has been calculated
    if(!("passedVariability" %in% colnames(mcols(fds, type="j")))){
        stop("Please calculate the expression filter values first with the ",
                "filterVariability function.")
    }
    
    # get plotting data
    dt <- data.table(
        value=pmax(mcols(fds, type="j")[['maxDPsi3']], 
                    mcols(fds, type="j")[['maxDPsi5']],
                    mcols(fds, type="j")[['maxDThetaDonor']],
                    mcols(fds, type="j")[['maxDThetaAcceptor']]),
        passed=mcols(fds, type="j")[['passedVariability']])
    if(isTRUE(onlyExpressedIntrons)){
        dt[,passed:=mcols(fds, type="j")[['passed']]]
    }
    
    # check if file with removed counts exists and add them when it exists
    nonVarDir <- file.path(workingDir(fds), "savedObjects", nameNoSpace(fds),
                            "nonVariableJunctions")

    if(dir.exists(nonVarDir)){
        nV_stored <- loadHDF5SummarizedExperiment(dir=nonVarDir) 
        nonVar_dt <- data.table(
            value=pmax(mcols(nV_stored)[['maxDPsi3']], 
                        mcols(nV_stored)[['maxDPsi5']],
                        mcols(nV_stored)[['maxDThetaDonor']],
                        mcols(nV_stored)[['maxDThetaAcceptor']]),
            passed=FALSE)
        dt <- rbind(dt, nonVar_dt)
    }
    
    dt[,passed:=factor(passed, levels=c(TRUE, FALSE))]
    colors <- brewer.pal(3, "Dark2")[seq_len(2)]
    ggplot(dt, aes(value, fill=passed)) +
        geom_histogram(bins=bins) +
        scale_y_log10() + 
        scale_fill_manual(values=colors, name="Passed",
                            labels=c("True", "False")) +
        xlab(bquote("Maximal Junction" ~ Delta*Psi[5] ~ "or" ~ Delta*Psi[3])) +
        ylab("Count") +
        ggtitle("Variability filtering") +
        theme_bw() +
        theme(legend.position=legend.position)
}


plotCountCorHeatmap.FRASER <- function(object,
                    type=psiTypes_avail, logit=FALSE, 
                    topN=50000, topJ=5000, minMedian=1, minCount=10, 
                    main=NULL, normalized=FALSE, show_rownames=FALSE,
                    show_colnames=FALSE, minDeltaPsi=0.1, annotation_col=NA,
                    annotation_row=NA, border_color=NA, nClust=5,
                    plotType=c("sampleCorrelation", "junctionSample"),
                    sampleClustering=NULL, plotMeanPsi=TRUE, plotCov=TRUE, ...){

    type <- match.arg(type)
    plotType <- match.arg(plotType)

    # use counts as matrix, otherwise x(fds,...) does not work later on
    counts(object, type=type, side="other", HDF5=FALSE)      <-
        as.matrix(counts(object, type=type, side="other"))
    counts(object, type=type, side="ofInterest", HDF5=FALSE) <-
        as.matrix(counts(object, type=type, side="ofInterest"))

    kmat <- K(object, type=type)
    nmat <- N(object, type=type)

    expRowsMedian <- rowMedians(kmat) >= minMedian
    expRowsMax    <- rowMax(kmat)     >= minCount
    table(expRowsMax & expRowsMedian)

    skmat <- kmat[expRowsMax & expRowsMedian,]
    snmat <- nmat[expRowsMax & expRowsMedian,]
    
    # check that we have at least 1 read for each sample
    minCovColSums <- colSums(snmat > minCount)
    if(any(minCovColSums < 2)){
        message("Warning:",
                " The following samples do not have at least 2 junctions",
                " with the minimum read coverage after filtering!",
                " They will be disregarded from the plot. ", 
                "\nAffected IDs are: \n\t", paste(collapse=", ", 
                        names(minCovColSums)[minCovColSums < 2]))
        ids2plot <- logical(length(minCovColSums))
        ids2plot[minCovColSums >= 2] <- TRUE
        skmat <- skmat[,ids2plot]
        snmat <- snmat[,ids2plot]
        object <- object[,ids2plot]
    }

    xmat <- (skmat + 1)/(snmat + 2)
    if(isTRUE(logit)){
        xmat <- qlogisWithCap(xmat)
    }
    xmat[snmat < minCount] <- NA
    xmat_rc    <- xmat - rowMeans(xmat, na.rm=TRUE)
    
    xmat_rc_sd <- rowSds(xmat_rc, na.rm=TRUE)
    nrNonNA <- rowSums(!is.na(xmat_rc))
    xmat_rc_sd[nrNonNA > 0.5*ncol(xmat_rc)] <- min(xmat_rc_sd)
    plotIdx <- rank(xmat_rc_sd) >= length(xmat_rc_sd) - topN
    xmat_rc_2_plot <- xmat_rc[plotIdx,]
    cormatS <- cor(xmat_rc_2_plot, use="pairwise", method="spearman")
    if(isTRUE(normalized)){
        pred_mu <- as.matrix(predictedMeans(object, type=type)[
            expRowsMax & expRowsMedian,][plotIdx,])
        if(isTRUE(logit)){
            pred_mu <- qlogisWithCap(pred_mu)
        }
        pred_mu[(snmat < minCount)[plotIdx,]] <- NA
        lpred_mu_rc <- pred_mu - rowMeans(pred_mu, na.rm=TRUE)
        xmat_rc_2_plot <- xmat_rc_2_plot - lpred_mu_rc
    }
    cormat <- cor(xmat_rc_2_plot, use="pairwise", method="spearman")
    breaks <- seq(-1, 1, length.out=50)

    if(plotType == "junctionSample"){

        if(isTRUE(normalized)){
            pred_mu <- as.matrix(predictedMeans(object, type=type)[
                expRowsMax & expRowsMedian,])
            if(isTRUE(logit)){
                pred_mu <- qlogisWithCap(pred_mu)
                breaks <- seq(-5, 5, length.out=50)
            }
            lpred_mu_rc <- pred_mu - rowMeans(pred_mu)
            xmat_rc <- xmat_rc - lpred_mu_rc
        }

        object <- object[expRowsMax & expRowsMedian,,by=type]
        j2keepVa <- variableJunctions(object, type, minDeltaPsi)
        j2keepDP <- rowQuantiles(kmat[expRowsMax & expRowsMedian,],
                                    probs=0.75) >= 10
        j2keep <- j2keepDP & j2keepVa
        xmat_rc_2_plot <- xmat_rc[j2keep,]
        mostVarKeep <- subsetKMostVariableJunctions(
                object[j2keep,,by=type], type, topJ)
        xmat_rc_2_plot <- xmat_rc_2_plot[mostVarKeep,]
        rownames(xmat_rc_2_plot) <- seq_len(nrow(xmat_rc_2_plot))

    }


    if(is.character(annotation_col)){
        annotation_col <- getColDataAsDFFactors(object, annotation_col)
    }
    if(is.character(annotation_row)){
        annotation_row <- getColDataAsDFFactors(object, annotation_row)
    }

    # annotate with sample clusters
    if(is.null(sampleClustering)){
        # annotate samples with clusters from sample correlation heatmap
        nClust <- min(nClust, nrow(cormatS))
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
            rowMeans(plogis(xmat), na.rm=TRUE)
        } else{
            rowMeans(xmat, na.rm=TRUE)
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
        main <- ifelse(normalized, "Normalized intron-centered ", 
                        "Raw intron-centered ")
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
#' Plot count correlation
#'
#' Count correlation heatmap function
#'
#' @rdname plotFunctions
#' @export
setMethod("plotCountCorHeatmap", signature="FraserDataSet", 
        plotCountCorHeatmap.FRASER)

#'
#' Plot coverage from bam files for given genomic range and sample ids
#'
#' @rdname plotFunctions
#' @export
plotBamCoverage <- function(fds, gr, sampleID,
        control_samples=sample(
            samples(fds[, which(samples(fds) != sampleID)]), 
            min(3, ncol(fds)-length(sampleID))), 
        txdb=NULL, min_junction_count=20, 
        highlight_range=NULL, highlight_range_color="firebrick", 
        color_annotated="gray", color_novel="goldenrod3", 
        color_sample_interest="firebrick", color_control_samples="dodgerblue4",
        toscale=c("exon", "gene", "none"), mar=c(2, 10, 0.1, 5),
        curvature_splicegraph=1, curvature_coverage=1, cex=1,
        splicegraph_labels=c("genomic_range", "id", "name", "none"),
        splicegraph_position=c("top", "bottom"), ...){
    
    if(missing(fds)){
        stop("Missing input: fds (FraserDataSet object)")
    } else{
        stopifnot(is(fds, "FraserDataSet"))
    }
    if(missing(gr)){
        stop("Missing input gr (genomic range to plot).")
    } else{
        stopifnot(is(gr, "GenomicRanges"))
        stopifnot(length(gr) > 0)
    }
    if(missing(sampleID)){
        stop("Missing input: sample_of_interest")
    }
    toscale <- match.arg(toscale)
    splicegraph_labels <- match.arg(splicegraph_labels)
    splicegraph_position <- match.arg(splicegraph_position)
    
    # extract bam info for sample ids to plot
    all_sids <- c(sampleID, control_samples)
    si_out <- getSGSeqSI(fds, all_sids)
    sgseq_si <- si_out[[1]]
    fds <- si_out[[2]]
    
    # collapse input ranges if several 
    if(any(strand(gr) == "*")){
        # seems to throw an error with * strand so guessing strand instead
        if(all(strand(gr) == "*")){
            guessStrand <- "+"
        } else{
            guessStrand <- strand(gr[strand(gr) != "*"])[1]
        }
        strand(gr)[strand(gr) == "*"] <- guessStrand
        warning("Input genomic ranges contained unstranded ranges.\n", 
                "This function needs strand information, guessing strand to ", 
                "be ", guessStrand, ".")
    }
    if(!all(strand(gr) == strand(gr[1,]))){
        warning("Input genomic ranges contained ranges on different strands,\n", 
                "only showing coverage for the ", strand(gr[1,]), " strand.")
        strand(gr) <- rep(strand(gr[1,]), length(gr))
    }
    gr <- range(gr)
    gr <- keepSeqlevels(gr, unique(as.character(seqnames(gr))))
    
    # convert highlight_range to GRangesList if not
    if(!is.null(highlight_range) && !is(highlight_range, "GRangesList")){
        stopifnot(is(highlight_range, "GRanges"))
        highlight_range <- GRangesList(highlight_range)
    }
    
    # extract splice graph
    sgfc_pred <- SGSeq::analyzeFeatures(sgseq_si, which = gr, 
                                min_junction_count=min_junction_count, psi=0, 
                                ...)
    
    # overlap detected junctions with annotation
    if(!is.null(txdb)){
        # subset to chr of interest
        seqlevels(txdb) <- unique(as.character(seqnames(gr)))
        
        # extract transcript features with SGSeq package
        txf <- SGSeq::convertToTxFeatures(txdb)
        txf <- txf[txf %over% gr]
        
        # restore seqlevels of txdb object
        seqlevels(txdb) <- seqlevels0(txdb)
        
        # annotate splice junctions with annotation features
        sgfc_pred <- SGSeq::annotate(sgfc_pred, txf)
    } else{
        # when no annotation is given, show everything in the same color
        color_novel <- color_annotated
    }
    
    # get genomic positions for first and last exon in given range
    if(splicegraph_labels == "genomic_range"){
        # tell plotSpliceGraph function to use custom labels
        splicegraph_labels <- "label"
        # create custom labels (only for first and last exon for readability)
        mcols(sgfc_pred)$label <- ""
        exons <- which(SGSeq::type(sgfc_pred) == "E" & 
                        rowRanges(sgfc_pred) %over% gr)
        exons <- unique(c(exons[1], tail(exons, n=1)))
        if(length(exons) == 1){
            mcols(sgfc_pred)$label[exons] <- 
                paste(seqnames(sgfc_pred), 
                        paste(start(sgfc_pred), end(sgfc_pred), sep="-"),
                        strand(sgfc_pred), sep=":")[exons]
        }
        if(length(exons) == 2){
            mcols(sgfc_pred)$label[exons[1]] <- 
                paste(seqnames(sgfc_pred), 
                        start(sgfc_pred), 
                        strand(sgfc_pred), sep=":")[exons[1]]
            mcols(sgfc_pred)$label[exons[2]] <- 
                paste(seqnames(sgfc_pred), 
                        end(sgfc_pred), 
                        strand(sgfc_pred), sep=":")[exons[2]]
        }
    }
    
    # plot splice graph and coverage of junctions from bam
    nr_sa2p <- length(all_sids)
    par(mfrow = c(nr_sa2p+1, 1), mar=mar, cex=cex) 
    if(splicegraph_position == "top"){
        SGSeq::plotSpliceGraph(rowRanges(sgfc_pred), 
                which=gr, 
                toscale=toscale, 
                color=color_annotated,
                color_novel=color_novel,
                ypos=c(0.25, 0.1),
                ranges=highlight_range,
                ranges_color=highlight_range_color, 
                ranges_ypos=c(0.01, 0.02),
                curvature=curvature_splicegraph,
                label=splicegraph_labels)
    }
    for (j in seq_along(sampleID)) {
        SGSeq::plotCoverage(
                sgfc_pred[, which(colnames(sgfc_pred) == sampleID[j])], 
                which = gr,
                toscale = toscale, 
                label=sampleID[j],
                color=color_sample_interest,
                curvature=curvature_coverage)
    }
    for (j in seq_along(control_samples)) {
        SGSeq::plotCoverage(
                sgfc_pred[, which(colnames(sgfc_pred) == control_samples[j])],
                which = gr,
                toscale = toscale, 
                label=control_samples[j],
                color=color_control_samples,
                curvature=curvature_coverage)
    }
    if(splicegraph_position == "bottom"){
        SGSeq::plotSpliceGraph(rowRanges(sgfc_pred), 
                which=gr, 
                toscale=toscale, 
                color_novel=color_novel,
                ypos=c(0.25, 0.1),
                ranges=highlight_range,
                ranges_color=highlight_range_color, 
                ranges_ypos=c(0.01, 0.02),
                curvature=curvature_splicegraph,
                label=splicegraph_labels)
    }
    
    return(invisible(fds))
}

#'
#' Plot coverage from bam files for given row of results table
#'
#' @rdname plotFunctions
#' @export
plotBamCoverageFromResultTable <- function(fds, result, show_full_gene=FALSE, 
                txdb=NULL, orgDb=NULL, res_gene_col="hgncSymbol",
                res_geneid_type="SYMBOL", txdb_geneid_type="ENTREZID",
                left_extension=1000, right_extension=1000, ...){
    stopifnot(is(fds, "FraserDataSet"))
    
    if(is(result, "GenomicRanges")){
        result <- as.data.table(result)
    }
    
    stopifnot(is.data.table(result))
    stopifnot(result[,.N] == 1)
    
    sid <- result[,sampleID]
    jidx <- getIndexFromResultTable(fds, result)
    outlier_range <- rowRanges(fds, type=result[,type])[jidx,]
    
    # showing either full range of the gene in which the outlier occured
    if(show_full_gene == TRUE){
        if(missing(txdb)){
            stop("Missing input: txdb (for extracting gene range)")
        }
        if(missing(orgDb)){
            stop("Missing input: orgDb (for mapping of IDs to txdb)")
        }
        result_gene <- result[,get(res_gene_col)]
        result_gene <- strsplit(result_gene, ";", fixed=TRUE)[[1]]
        if(is.data.table(orgDb)){
            tmp <- merge(x=as.data.table(genes(txdb))[,.(gene_id)], y=orgDb, 
                         by.y=txdb_geneid_type, by.x="gene_id", all.x=TRUE, 
                         sort=FALSE)[,.(gene_id, feature=get(res_geneid_type))]
            setnames(tmp, "feature", res_geneid_type)
            txdb_geneid <- tmp[get(res_geneid_type) %in% result_gene, gene_id]
        } else {
            tmp <- as.data.table(
                select(orgDb, 
                       keys=result_gene, 
                       columns=txdb_geneid_type, 
                       keytype=res_geneid_type)
            )
            txdb_geneid <- tmp[, get(txdb_geneid_type)]
        }
        gr <- genes(txdb, filter=list("gene_id"=txdb_geneid))
        if(length(gr) == 0){
            stop("Could not extract genomic coordinates for input gene.")
        }
    } else{
        # or just showing a certain region around the outlier junction
        gr <- outlier_range
        start(gr) <- start(gr) - left_extension
        end(gr) <- end(gr) + right_extension
    }
    
    # if several genes overlap, only show those on same strand as outlier
    if(as.character(strand(outlier_range)) != "*" & 
        length(gr[strand(gr) == strand(outlier_range),]) > 0){
        gr <- gr[strand(gr) == strand(outlier_range),]
    }
    
    # create the coverage plot for the given outlier
    fds <- plotBamCoverage(fds, 
                            gr=gr, 
                            sampleID=sid, 
                            txdb=txdb,
                            highlight_range=outlier_range, 
                            ...)
    return(invisible(fds))
}

plotManhattan.FRASER <- function(object, sampleID, 
                                type=psiTypes_avail, 
                                main=paste0("sampleID = ", sampleID), 
                                color_chr=c("black", "darkgrey"),
                                ...){
    # check arguments
    stopifnot(is(object, "FraserDataSet"))
    stopifnot(sampleID %in% samples(object))
    type <- match.arg(type)
    additional_args <- list(...)
    padjCutoff <- 0.05
    if("padjCutoff" %in% names(additional_args)){
        padjCutoff <- additional_args$padjCutoff
    }
    deltaPsiCutoff <- 0.3
    if("deltaPsiCutoff" %in% names(additional_args)){
        deltaPsiCutoff <- additional_args$deltaPsiCutoff
    }
    
    # extract neccessary informations
    gr_sample <- rowRanges(object, type=type)
    seqlevelsStyle(gr_sample) <- seqlevelsStyle(object)
    mcols(gr_sample)[,"pvalue"] <- -log10(
        pVals(object, type=type, level="junction")[,sampleID])
    mcols(gr_sample)[,"padjust"] <- -log10(
        padjVals(object, type=type, level="site")[,sampleID])
    mcols(gr_sample)[,"delta"] <- deltaPsiValue(object, type=type)[,sampleID]
    
    # only one point per donor/acceptor site (relevant only for psi5 and psi3)
    index <- FRASER:::getSiteIndex(object, type=type)
    nonDup <- !duplicated(index)
    gr_sample <- gr_sample[nonDup,]
    
    # Sort granges for plot
    gr_sample <- sortSeqlevels(gr_sample)
    gr_sample <- sort(gr_sample)
    
    # find outlier indices
    if(!type %in% c("psi3", "psi5")){
        outlier_idx <- which(gr_sample$padjust >= -log10(padjCutoff) &
                                abs(gr_sample$delta) >= deltaPsiCutoff)
    } else{
        outlier_idx <- which(gr_sample$padjust >= -log10(padjCutoff))
    }
    message("highlighting ", length(gr_sample[outlier_idx,]), " outliers ...")
    
    # plot manhattan plot
    plotGrandLinear.adapted(gr_sample, aes(y=pvalue), 
                            color=color_chr, 
                            highlight.gr=gr_sample[outlier_idx,],
                            highlight.overlap="equal") + 
        labs(title=main)
    
}

#'
#' Plot manhattan plot of junction pvalues
#'
#' @rdname plotFunctions
#' @export
setMethod("plotManhattan", signature="FraserDataSet", 
          plotManhattan.FRASER)


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

#' 
#' used to cap the qlogis for the correlation heatmap
#' 
#' @noRd
qlogisWithCap <- function(x, digits=2){
    x <- round(x, digits)
    x <- pmin(pmax(x, 10^-digits), 1-10^-digits)
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
ggplotLabelPsi <- function(type, asCharacter=FALSE){
    type <- as.character(type)
    if(isFALSE(asCharacter)){
        vapply(type, FUN=function(x)
            switch (x,
                    jaccard = c(bquote(jaccard~intron~index)),
                    psi5 = c(bquote(psi[5])),
                    psi3 = c(bquote(psi[3])),
                    theta = c(bquote(theta))),
            FUN.VALUE=c(bquote(psi[3])))
    } else{
        vapply(type, FUN=function(x)
            switch (x,
                    jaccard = "Intron-Jaccard-Index",
                    psi5 = "psi[5]",
                    psi3 = "psi[3]",
                    theta = "theta"),
            FUN.VALUE=character(1))
    }
}

#'
#' Extract info from bam files needed for SGSeq functions to work
#' 
#' @noRd
getSGSeqSI <- function(fds, sample_ids){
    
    # check if bam info is already stored in fds for given samples
    if("SGSeq_sampleinfo" %in% names(metadata(fds))){
        si <- metadata(fds)[["SGSeq_sampleinfo"]]
        si <- si[si$sample_name %in% sample_ids,]
        if(nrow(si) != length(sample_ids)){
            # add bam info for missing sample_ids
            missing_ids <- sample_ids[!sample_ids %in% si$sample_name]
            message("Extracting SGSeq sample info from BAM files for samples ",
                    paste(missing_ids, collapse=", "), " ...")
            df_missing <- data.frame(
                sample_name=samples(fds)[samples(fds) %in% missing_ids],
                file_bam=bamFile(fds)[samples(fds) %in% missing_ids])
            si_new <- SGSeq::getBamInfo(df_missing, yieldSize=1e6)
            si_new$lib_size <- 50e6 # dummy value to speed up this part
            si <- rbind(si, si_new)
            metadata(fds)[["SGSeq_sampleinfo"]] <- 
                rbind(metadata(fds)[["SGSeq_sampleinfo"]], si_new)
        }
        return(list(si, fds))
    } else{
        message("Extracting SGSeq sample info from BAM files for samples ",
                paste(sample_ids, collapse=", "), " ...")
        df <- data.frame(
            sample_name=samples(fds)[samples(fds) %in% sample_ids],
            file_bam=bamFile(fds)[samples(fds) %in% sample_ids])
        si <- SGSeq::getBamInfo(df, yieldSize=1e6)  
        si$lib_size <- 50e6 # dummy value to speed up this part
        metadata(fds)[["SGSeq_sampleinfo"]] <- si
        return(list(si, fds))
    }
}

#'
#' Adapted function from ggbio package to create manhattan plot. 
#' Adapted to allow highlighting only ranges that exactly match. Uses functions 
#' from package biovizBase.
#'
#' @noRd
plotGrandLinear.adapted <- function (obj, ..., facets, space.skip = 0.01, 
        geom = NULL, cutoff = NULL, cutoff.color = "red", cutoff.size = 1, 
        legend = FALSE, xlim, ylim, xlab, ylab, main, highlight.gr = NULL, 
        highlight.name = NULL, highlight.col = "red", highlight.label = TRUE, 
        highlight.label.size = 5, highlight.label.offset = 0.05, 
        highlight.label.col = "black", 
        highlight.overlap = c("any", "start", "end", "within", "equal"),
        spaceline = FALSE){
    if (is.null(geom)) 
        geom <- "point"
    args <- list(...)
    args.aes <- biovizBase::parseArgsForAes(args)
    args.non <- biovizBase::parseArgsForNonAes(args)
    two.color <- c("#0080FF", "#4CC4FF")
    .is.seq <- FALSE
    if (!"colour" %in% names(args.aes)) {
        if (!any(c("color", "colour") %in% names(args.non))) {
            .color <- two.color
            args.aes$color <- as.name("seqnames")
            .is.seq <- TRUE
        }
        else {
            if (length(args.non$color) > 1) {
                .color <- args.non$color
                args.aes$color <- as.name("seqnames")
                .is.seq <- TRUE
                args.non <- args.non[!names(args.non) %in% c("colour", 
                                                             "color")]
            }
        }
    }
    else {
        if (quo_name(args.aes$colour) == "seqnames") 
            args.aes$colour <- as.name("seqnames")
    }
    if (!"y" %in% names(args.aes)) 
        stop("need to provide y")
    args.non$coord <- "genome"
    args.non$space.skip <- space.skip
    args.non$geom <- geom
    args.non$object <- obj
    aes.res <- do.call(aes, args.aes)
    p <- do.call(ggbio::autoplot, c(list(aes.res), args.non))
    if (!legend) 
        p <- p + theme(legend.position = "none")
    if (!missing(ylab)) 
        p <- p + ylab(ylab)
    if (!is.null(cutoff)) 
        p <- p + geom_hline(yintercept = cutoff, color = cutoff.color, 
                            size = cutoff.size)
    chrs <- names(seqlengths(obj))
    if (.is.seq) {
        N <- length(chrs)
        cols <- rep(.color, round(N/length(.color)) + 1)[1:N]
        names(cols) <- chrs
        p <- p + scale_color_manual(values = cols)
    }
    if (!missing(facets)) {
        args$facets <- facets
        args.facets <- subsetArgsByFormals(args, facet_grid, 
                                           facet_wrap)
        facet <- .buildFacetsFromArgs(obj, args.facets)
        p <- p + facet
    }
    p <- p + theme(panel.grid.minor = element_blank())
    if (!is.null(highlight.gr)) {
        highlight.overlap <- match.arg(highlight.overlap)
        idx <- findOverlaps(obj, highlight.gr, type=highlight.overlap)
        .h.pos <- lapply(split(queryHits(idx), subjectHits(idx)), 
                         function(id) {
                             gr <- GRanges(as.character(seqnames(p@data))[id][1], 
                                           IRanges(start = min(start(p@data[id])), end = max(end(p@data[id]))))
                             val <- max(as.numeric(values(p@data[id])[, quo_name(args.aes$y)]))
                             val <- val * (1 + highlight.label.offset)
                             values(gr)$val <- val
                             gr
                         })
        .h.pos <- suppressWarnings(do.call("c", unname(.h.pos)))
        if (length(.h.pos)) {
            if (is.null(highlight.name)) {
                highlight.name <- names(highlight.gr)
            }
            else {
                highlight.name <- values(highlight.gr)[, highlight.name]
            }
            p <- p + geom_point(data = mold(p@data[queryHits(idx)]), 
                                do.call(aes, list(x = substitute(midpoint), y = args.aes$y)), 
                                color = highlight.col)
            if (!is.null(highlight.name)) {
                seqlevels(.h.pos, pruning.mode = "coarse") <- seqlevels(obj)
                suppressWarnings(seqinfo(.h.pos) <- seqinfo(obj))
                .trans <- transformToGenome(.h.pos, space.skip = space.skip)
                values(.trans)$mean <- (start(.trans) + end(.trans))/2
                values(.trans)$names <- highlight.name
                p <- p + geom_text(data = mold(.trans), size = highlight.label.size, 
                                   vjust = 0, color = highlight.label.col, do.call(aes, 
                                                                                   list(x = substitute(mean), y = as.name("val"), 
                                                                                        label = as.name("names"))))
            }
        }
    }
    if (spaceline) {
        vline.df <- p@ggplot$data
        vline.df <- do.call(rbind, by(vline.df, vline.df$seqnames, 
                                      function(dd) {
                                          data.frame(start = min(dd$start), end = max(dd$end))
                                      }))
        gap <- (vline.df$start[-1] + vline.df$end[-nrow(vline.df)])/2
        p <- p + geom_vline(xintercept = gap, alpha = 0.5, color = "gray70") + 
            theme(panel.grid = element_blank())
    }
    if (!missing(main)) 
        p <- p + labs(title = main)
    if (!missing(xlim)) 
        p <- p + xlim(xlim)
    if (!missing(ylim)) 
        p <- p + ylim(ylim)
    if (missing(xlab)) 
        xlab <- ""
    p <- p + ggplot2::xlab(xlab)
    p
}
