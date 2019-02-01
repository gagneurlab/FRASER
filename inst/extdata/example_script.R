devtools::document()
devtools::load_all()
devtools::install(".")
unloadNamespace("FraseR")

library("FraseR")

bpparam <- SnowParam(workers=8, progressbar=TRUE)
bpparam <- SerialParam()
register(bpparam)

fds <- loadFraseRDataSet("~/projects/data", "kremer-bader-et-al", upgrade=TRUE)
probE <- 0.01
featureExclusionMask(fds) <- sample(c(TRUE, FALSE), nrow(fds), replace=TRUE, prob=c(probE,1-probE))
table(featureExclusionMask(fds))
debug(fitAutoencoder)
fds <- fitAutoencoder(fds, q=24, verbose=TRUE, BPPARAM=bpparam)


