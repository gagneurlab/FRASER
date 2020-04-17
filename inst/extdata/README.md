# Datasets accompanying the FRASER package

The FRASER packages contains an example dataset that is a subset of the 
[Kremer et al study](https://doi.org/10.1038/ncomms15824). The main 
purpose of those files are to demonstrate the functionality of the 
FRASER package.

## BAM files

The three BAM files in `inst/extdata/bam` belong to one of the diagnost
individuals from the study. In short, they contain Illumina RNA-seq data
from fibroblast samples and aligned with STAR against the UCSC hg19 assembly.
For this package only reads spanning the genes `MCOLN1` and `TIMMDC1` are
extracted. More details on the protocol can be found in the 
[Method section](https://doi.org/10.1038/ncomms15824).

## Count matrices 

The split-read and splice-site-read count matrices
(`inst/extdata/raw_junction_counts.tsv.gz` and 
`inst/extdata/raw_site_counts.tsv.gz`) are also a subset of the 
[Kremer et al study](https://doi.org/10.1038/ncomms15824). Here, the
individuals with diagnosed splice defects and some controls are picked.
The counts are extracted from the BAM files with the function 
`FRASER::countRNAData`. Again it contains only a subset of count for 
the genes `TIMMDC1` and `MCOLN1`.

## Sample annotation

Both datasets are accompanied by a pseudo anonymised sample annotation. 

