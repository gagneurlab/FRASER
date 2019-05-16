annotate_strand_by_seq <- function(fds, genome, ...){
    #load genome
    if (genome=="hg19") {
        library(BSgenome.Hsapiens.UCSC.hg19)
        genome <- BSgenome.Hsapiens.UCSC.hg19
    }
    else if (genome=="GRCh37"){
        library(AnnotationHub)
        genome <- ah[["AH134"]]
    }
    else if (genome=="GRCh38"){
        library(BSgenome.Hsapiens.NCBI.GRCh38)
        genome <- BSgenome.Hsapiens.NCBI.GRCh38
    }


    # get ranges of interest
    gr <- granges(fds)

    # get the dinucleotide seq
    gr_start <- resize(gr, 2, fix = "start")
    gr_end <- resize(gr, 2, fix = "end")

    mcols(fds, type="j")[["psi5_dinuc"]] <- getSeq(genome, gr_start, as.character = TRUE)
    mcols(fds, type="j")[["psi3_dinuc"]] <- getSeq(genome, gr_end, as.character = TRUE)

    # set the new seq based annotation
    pos <- mcols(fds, type="j")[["psi5_dinuc"]] %in% c("GT", "GC") & mcols(fds, type="j")[["psi3_dinuc"]] == "AG"
    neg <- mcols(fds, type="j")[["psi5_dinuc"]] == "CT" & mcols(fds, type="j")[["psi3_dinuc"]] %in% c("AC", "GC")


    mcols(fds, type="j")[["predicted_strand"]] <- "*"
    mcols(fds, type="j")[pos, ][["predicted_strand"]] <- "+"
    mcols(fds, type="j")[neg, ][["predicted_strand"]] <- "-"

    # return object
    return(fds)
}


strand_summary <- function(fds, ...){
    GTpos <- mcols(fds, type="j")[["psi5_dinuc"]] == "GT" & mcols(fds, type="j")[["psi3_dinuc"]] == "AG"
    GCpos <- mcols(fds, type="j")[["psi5_dinuc"]] == "GC" & mcols(fds, type="j")[["psi3_dinuc"]] == "AG"

    GTneg <- mcols(fds, type="j")[["psi5_dinuc"]] == "CT" & mcols(fds, type="j")[["psi3_dinuc"]] == "AC"
    GCneg <- mcols(fds, type="j")[["psi5_dinuc"]] == "CT" & mcols(fds, type="j")[["psi3_dinuc"]] == "GC"

    pos <- mcols(fds, type="j")[["psi5_dinuc"]] == "GT" | mcols(fds, type="j")[["psi3_dinuc"]] == "AG"
    neg <- mcols(fds, type="j")[["psi5_dinuc"]] == "CT" | mcols(fds, type="j")[["psi3_dinuc"]] == "AC"
    posneg <- pos&neg

    half_pos <- xor(mcols(fds, type="j")[["psi5_dinuc"]] == "GT", mcols(fds, type="j")[["psi3_dinuc"]] == "AG") & !GCpos & !posneg
    half_neg <- xor(mcols(fds, type="j")[["psi5_dinuc"]] == "CT", mcols(fds, type="j")[["psi3_dinuc"]] == "AC") & !GCneg & !posneg

    none <- !(pos|neg)

    strand_stat <- c(GTpos = sum(GTpos), GCpos = sum(GCpos), GTneg = sum(GTneg), GCneg = sum(GCneg),
                     half_pos = sum(half_pos), half_neg = sum(half_neg), posneg = sum(posneg), none = sum(none), total = length(fds)
    )
    return(strand_stat)
}


