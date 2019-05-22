devtools::document()
devtools::load_all()

bpparam <- MulticoreParam(50, 150, progressbar=TRUE)
bpparam <- MulticoreParam(30, 90, progressbar=TRUE)
bpparam <- SnowParam(workers=8, progressbar=TRUE)
bpparam <- SerialParam()
register(bpparam)

# real data
fdsf <- loadFraseRDataSet("/dev/shm/fraser-mertes", "kremer-bader-et-al")
fds <- fdsf
probE <- 0.1
q <- 20

# real data
fds <- loadFraseRDataSet("/s/project/fraser/analysis/datasets", "prokisch_batch5")
fds <- loadFraseRDataSet("~/projects/data", "kremer-bader-et-al", upgrade=TRUE)
probE <- 0.1
q <- 24

# toy data
if(FALSE){
    fdse <- makeExampleFraseRDataSet()
    fdse <- saveFraseRDataSet(fdse)
}
fdse <-loadFraseRDataSet("~/FraseR", "Data_Analysis")
fds <- fdse
probE <- 0.5
q <- 3

debug(fitAutoencoder)
debug(subsetFraseR)
debug(singleDFit)
debug(updateD)

type <- "psiSite"
for(type in c("psi5", "psi3", "psiSite")){
    # set current type
    currentType(fds) <- type

    # subset fitting
    featureExclusionMask(fds) <- sample(c(TRUE, FALSE),
            nrow(mcols(fds, type=type)), replace=TRUE, prob=c(probE,1-probE))
    print(table(featureExclusionMask(fds)))

    # run autoencoder
    fds <- fitAutoencoder(fds, q=q, type=type, verbose=TRUE, BPPARAM=bpparam, iterations=15)

    # calculate stats
    fds <- calculateZscore(fds, type=type)
    fds <- calculatePvalues(fds, type=type)
    fds <- calculatePadjValues(fds, type=type)
}

fds <- saveFraseRDataSet(fds)

fds <- annotateRanges(fds)
results <- results(fds, fdrCut=0.9)


granges(rowRanges(fds))

# TIMMDC1
grList <- unlist(GRangesList(
    TIMMDC1_3C = GRanges("chr3", IRanges(119234787, 119236051), type="psi3", samples="MUC1344|MUC1365"),
    TIMMDC1_3A = GRanges("chr3", IRanges(119232567, 119236051), type="psi3", samples="MUC1344|MUC1365"),
    TIMMDC1_5C = GRanges("chr3", IRanges(119234787, 119236051), type="psi5", samples="MUC1344|MUC1365"),
    TIMMDC1_5A = GRanges("chr3", IRanges(119232567, 119234706), type="psi5", samples="MUC1344|MUC1365"),
    MCOLN1_S1  = GRanges("chr19", IRanges(7592514, 7592515), type="psiSite", samples="MUC1361"),
    MCOLN1_S2  = GRanges("chr19", IRanges(7592749, 7592750), type="psiSite", samples="MUC1361"),
    CLPP_5C    = GRanges("chr19", IRanges(6361463, 6368919), type="psi3", samples="MUC1350"),
    PTPRS	   = GRanges("chr19", IRanges(5286246, 5340674), type="psi5", samples="MUC1345"),
    SFXN4	   = GRanges("chr10", IRanges(120921851, 120921852), type="psiSite", samples="76624"),
    LPIN1	   = GRanges("chr2", IRanges(11955367, 11959609), type="psi5", samples="MUC1391"),
    PANK2	   = GRanges("chr20", IRanges(3888926, 3891223), type="psi3", samples="MUC1393"),
    TALDO1	   = GRanges("chr11", IRanges(760254, 763343), type="psi5", samples="MUC1410"),
    COASY	   = GRanges("chr17", IRanges(40714238, 40715978), type="psi3", samples="MUC1395"),
    DPYD	   = GRanges("chr1", IRanges(98198027, 98205947), type="psi3", samples="MUC1423"),
    TAZ	       = GRanges("chrX", IRanges(153641880, 153641881), type="psiSite", samples="MUC1398")
))
grList

gr <- grList[15]
ovtype <- "equal"
for(gr in grList){
    gr
    type <- gr$type
    currentType(fds) <- type
    sampleIDs <- unlist(strsplit(gr$samples, "\\|"))

    fdsgr <- granges(rowRanges(fds, type=type))

    idx <- to(findOverlaps(gr, fdsgr, type = ovtype))
    as.vector(pVals(fds)[idx,sampleIDs])
    sort(as.vector(pVals(fds)[idx,]))

    par(cex=2)
    plotJunction(idx[1], fds, dist="BetaBinomial", threshold = 0.1)
}

type <- "psiSite"
type <- "psi3"
pdf(paste0('kremer_ae_dpyd', type, '.pdf'), width=21, height=14)
#for(i in sample(1:nrow(fds), 50)){
for(i in which(mcols(fds, type="psi3")$hgnc_symbol == "DPYD")){
    currentType(fds) <- type
    plotJunction(i, fds, type=type)
}
dev.off()



getLine <- function(dt, idx){
    cat(dt[idx,paste0(hgnc_symbol, '\t= GRanges("', seqnames, '", IRanges(', start, ", ", end, '), type="', type, '", samples="', sampleID, '")')])
}

library(ggplot2)

# presentation plots
#
# # volcano TIMMDC1 / MUC1344
ggplot(as.data.table(results)[sampleID == "MUC1344"], aes(deltaPsi, -log10(abs(pvalue)))) + geom_point()
ggplot(as.data.table(results)[sampleID == "MUC1365"], aes(deltaPsi, -log10(abs(pvalue)))) + geom_point()


dt <- as.data.table(results)

filtdt <- dt[abs(p.adj) <= 0.1 &
            (expression + expressionOther >= 5) &
            meanTotalCts < 100 * expressionOther]

numEventsBySample <- filtdt[abs(deltaPsi) > 0.1, .N,by=c("sampleID", "hgnc_symbol")][
                        ,.N,by="sampleID"][,.(myrank=rank(N), N1=N, sampleID)][order(myrank)]
ggplot(numEventsBySample[,.(x=1:.N, y=N1)], aes(x=x, y=y)) +
    geom_histogram(stat="identity") +
    scale_y_log10() +
    ylab("Number of \naberrantly spliced genes") +
    xlab("Ranked subjects")


hist(filtdt[, deltaPsi], breaks=100)


# global qq plot
myPvals <- as.matrix(pVals(fds, type="psi5"))
myPvalsdt <- data.table(e=ppoints(prod(dim(myPvals))), o=sort(as.vector(abs(myPvals))))
myPvalsdt[o != 1]
ggplot(myPvalsdt[o < 0.7], aes(x=-log10(e), y=-log10(o))) +
    geom_hex() +
    geom_abline(slope=1, intercept=0, col="firebrick")



# volcano plot
dt2plot <- data.table(p=pVals(fds, type="psi5")[,"MUC1344"],
                      z=zScores(fds, type="psi5")[,"MUC1344"])
minPval <- min(-log10(abs(dt[sampleID == "MUC1344" & p.adj <= 0.1, (pvalue)])))
ggplot(dt2plot, aes(x=z, y=-log10(abs(p)))) +
    stat_binhex(aes(fill=log10(..count..))) +
    ylab("-log10(pvalue)") +
    xlab("Zscore") +
    geom_hline(yintercept=minPval, col="firebrick")
    #geom_vline(xintercept=c(-3, 3), col="firebrick")


dt2plot[order(abs(p))]


#
# TAZ
#
idx <- 368110
samEx <- N(fds)[idx,] != 0
plotJunction(idx, fds[,samEx], dist="BetaBinomial", threshold = 0.1)
