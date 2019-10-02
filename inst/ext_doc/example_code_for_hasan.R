#
# https://i12g-gagneurweb.informatik.tu-muenchen.de/gitlab/mertes/FraseR.git
#
#

sampleIDs <- unique(strsplit("
SRR1070159
SRR1070159
SRR1070184
SRR1070184
SRR1070208
SRR1070208
SRR1070232
SRR1070232
SRR1070260
SRR1070260
SRR1070286
SRR1070286
SRR1070308
SRR1070308
SRR1070332
SRR1070332
SRR1070358
SRR1070358
SRR1070382
SRR1070382
SRR1070403
SRR1070403
SRR1070429
SRR1070429
SRR1070455
SRR1070455
SRR1070479
SRR1070479
SRR1070503
SRR1070503
SRR1070525
SRR1070525
SRR1070549
SRR1070549
SRR1070573
SRR1070573
SRR1070597
SRR1070597
SRR1070620
SRR1070620
", "\n")[[1]][2:30])

# fds <- makeExampleFraseRDataSet()
# fds <- loadFraseRDataSet("~/FraseR/")

# load package
devtools::load_all("./")
library(data.table)

# create sample annotation and file path
anno <- data.table(sampleID=sampleIDs)
anno[,bamFile:=paste0("/s/project/sra-download/bamfiles/", sampleID, ".bam")]
anno


# create fraser object
fds <- FraseRDataSet(anno)
strandSpecific(fds) <- TRUE
fds

# set multicore params
parallel(fds) <- MulticoreParam(10, nrow(anno))

# setup the path where to store the data
workingDir(fds) <- "/s/project/fraser/analysis/testHasan"
name(fds) <- "GTEx1"

# count the data
fds <- countRNAData(fds, NcpuPerSample=4)

# compute PSI and N
fds <- calculatePSIValues(fds)


# save the object to retrieve it later again
fds <- saveFraseRDataSet(fds)



# filter if needed
# subsets directly the object
fds<- filterExpression(fds, minExpression=5, quantile=0.7, quantielMinExpression=0, filter=TRUE)


# retrive counts
type <- "psi5" # "psi3"

psiVals <- assay(fds, type)
k <- counts(fds, type=type, side="ofInterest")
n_minus_k <- counts(fds, type="psi5", side="otherSide")
n <- k + n_minus_k

dtout <- as.data.table(granges(rowRanges(fds, type=type)))
k
n
psiVals

# write out result file
dtkout <- cbind(dtout, as.data.table(k))
write.table(dtkout, row.names=FALSE, file="./tmpfileXXX.tsv")

