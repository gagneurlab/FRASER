context("FraseRDataSet methods")

tmp_settings <- createTestFraseRSettings()[3]

tmp_dataset  <- countRNAData(tmp_settings)

test_that("annotateRanges", {
    annotation_dataset <- annotateRanges(tmp_dataset)
    expect_that(mcols(annotation_dataset@splitReads)$hgnc_symbol, 
                equals(rep(c("TMEM39A", "TIMMDC1", "ARHGEF18", "PCP2", "MCOLN1", "PNPLA6"), c(1,16,2,2,15,2)))
    )
})



# system.time(a <- featureDT[1:1000,list(feature=paste(sort(unique(feature)), collapse=";")),by="from"])
# a
# featureDT[1:1000]
# 
# tmp_dt <- data.table(a=c(1:12), b=rep(1:4,3))
# unique(tmp_dt)
# tmp_dt <- setkey(tmp_dt, b)
# unique(tmp_dt, by="b")
# ?unique.data.table
# library(data.table)