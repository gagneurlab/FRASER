#'
#' get links
#'
#'
createFullLinkTable <- function(data, addHRef){
    if(is.null(data) | length(data) < 1)
        return(data.table())

    selSymbols <- data[,hgnc_symbol]
    selChr     <- data[,seqnames]
    selStart   <- data[,start]
    selEnd     <- data[,end]

    data.table(
        Symbol=selSymbols,
        Chromosome=selChr,
        Start=selStart,
        End=selEnd,
        GoogleSearch=getGoogleLink(selSymbols, addHRef),
        HGNC=getHGNCLink(selSymbols, addHRef),
        GeneCards=getGeneCardsLink(selSymbols, addHRef),
        Ensembl=getEnsemblSymbolLink(selSymbols, addHRef),
        OMIM=getOmimLink(selSymbols, addHRef),
        ExAC=getExACLink(selChr,selStart,selEnd, addHRef),
        genomAD=getGenomAD(selChr,selStart,selEnd, addHRef),
        EnsemblBrowser=getEnsemblBrowser(selChr,selStart,selEnd, addHRef)
    )
}

addHRef <- function(link, text, add=TRUE){
    if(!add){
        link
    } else {
        paste0('<a target="_blank" href="', link, '">', text, '</a>')
    }
}

getGoogleLink <- function(str, asHRef=TRUE){
    addHRef(paste0("https://www.google.com/search?q=", str), str, asHRef)
}

getHGNCLink <- function(symbol, asHRef=TRUE){
    addHRef(
        paste0(
            "http://www.genenames.org/cgi-bin/gene_symbol_report?match=",
            symbol
        ),
        symbol,
        asHRef
    )
}

getGeneCardsLink <- function(symbol, asHRef=TRUE) {
    addHRef(
        paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=", symbol),
        symbol,
        asHRef
    )
}

getEnsemblSymbolLink <- function(symbol, asHRef=TRUE){
    addHRef(getGoogleLink(paste("ENSEMBL", symbol), FALSE), symbol, asHRef)
}

getEnsemblBrowser <- function(chr, start, end, asHRef=TRUE, extend=1000){
    chr <- gsub("^chr", "", chr)
    addHRef(
        paste0("http://grch37.ensembl.org/Homo_sapiens/Location/View?r=",
            chr , "%3A", max(0, start-extend), "-", end+extend
        ),
        paste0(chr, ":", start, "-", end),
        asHRef
    )
}

getOmimLink <- function(symbol, asHRef=TRUE){
    addHRef(getGoogleLink(paste("OMIM", symbol), FALSE), symbol, asHRef)
}

getExACLink <- function(chr, start, end, asHRef=TRUE){
    chr <- gsub("^chr", "", chr)
    addHRef(
        paste0("http://exac.broadinstitute.org/region/",
            chr, "-", start, "-", end
        ),
        paste0(chr, ":", start, "-", end),
        asHRef
    )
}

getGenomAD <- function(chr, start, end, asHRef=TRUE){
    chr <- gsub("^chr", "", chr)
    addHRef(
        paste0("http://gnomad.broadinstitute.org/region/",
            chr, "-", start, "-", end
        ),
        paste0(chr, ":", start, "-", end),
        asHRef
    )
}
