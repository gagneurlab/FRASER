devtools::document()
devtools::load_all()

bpparam <- SnowParam(workers=8, progressbar=TRUE)
bpparam <- SerialParam()
register(bpparam)

# real data
fds <- loadFraseRDataSet("~/projects/data", "kremer-bader-et-al", upgrade=TRUE)
probE <- 0.1
q <- 24

# toy data
fds <- makeExampleFraseRDataSet()
probE <- 0.5
q <- 3

debug(fitAutoencoder)
debug(subsetFraseR)
debug(singleDFit)
debug(updateD)

type <- "psiSite"
for(type in c("psi5", "psi3", "psiSite")){

    # subset fitting
    featureExclusionMask(fds, type=type) <- sample(c(TRUE, FALSE),
            nrow(mcols(fds, type=type)), replace=TRUE, prob=c(probE,1-probE))
    table(featureExclusionMask(fds, type="psiSite"))

    # run autoencoder
    fds <- fitAutoencoder(fds, q=q, type=type, verbose=TRUE, BPPARAM=bpparam, iterations=5)

    # calculate stats
    fds <- calculateZscore(fds, type=type)
    fds <- calculatePvalues(fds, type=type)
}




