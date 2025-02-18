context("Test updateSeqlevelsStyle function")

test_that("hg38, UCSC to NCBI", {
  genome <- getBSgenome("hg38")
  bsgenome <- genome
  
  seqlevelsStyle(bsgenome) <- "NCBI"
  genome <- updateSeqlevelsStyle(genome, metadata(genome)$genome, "NCBI", metadata(genome)$provider)
  
  expect_equal(seqnames(genome), seqnames(bsgenome))
})

test_that("hg38, NCBI to UCSC", { 
  genome <- getBSgenome("GRCh38")
  bsgenome <- genome
  
  seqlevelsStyle(bsgenome) <- "UCSC"
  genome <- updateSeqlevelsStyle(genome, metadata(genome)$genome, "UCSC", metadata(genome)$provider)
  
  expect_equal(seqnames(genome), seqnames(bsgenome))
})

test_that("hg19, NCBI to UCSC", {
  genome <- getBSgenome("hs37d5")
  bsgenome <- genome
  
  seqlevelsStyle(bsgenome) <- "UCSC"
  genome <- updateSeqlevelsStyle(genome, metadata(genome)$genome, "UCSC", metadata(genome)$provider)
  
  expect_equal(seqnames(genome), seqnames(bsgenome))
})

test_that("hg19, UCSC to NCBI", {
  genome <- getBSgenome("hg19")
  bsgenome <- genome
  
  seqlevelsStyle(bsgenome) <- "NCBI"
  genome <- updateSeqlevelsStyle(genome, metadata(genome)$genome, "NCBI", metadata(genome)$provider)
  
  expect_equal(seqnames(genome), seqnames(bsgenome))
})
