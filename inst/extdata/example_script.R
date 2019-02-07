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
q <- 24

# real data
fds <- loadFraseRDataSet("~/projects/data", "kremer-bader-et-al", upgrade=TRUE)
probE <- 0.1
q <- 24

# toy data
fdse <- makeExampleFraseRDataSet()
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
    fds <- fitAutoencoder(fds, q=q, type=type, verbose=TRUE, BPPARAM=bpparam, iterations=10)

    # calculate stats
    fds <- calculateZscore(fds, type=type)
    fds <- calculatePvalues(fds, type=type)
    fds <- calculatePadjValues(fds, type=type)
}

rowRanges(fds)
granges(rowRanges(fds))

# TIMMDC1
library(GenomicRanges)
grList <- GRangesList(
    TIMMDC1 = GRanges("chr3", IRanges(119217368, 119243937), type="psi3"),
    MCOLN1  = GRanges("chr19", IRanges(7587496, 7598895), type="psiSite")
)
grList

gr <- grList[[1]]
for(gr in grList){
    type <- gr$type
    currentType(fds) <- type

    idx <- to(findOverlaps(gr, granges(rowRanges(fds, type=type))))
    hist(log10(as.vector(pVals(fds)[idx,c("MUC1344", "MUC1365")])))

    plotJunction(52984, fds)
}
