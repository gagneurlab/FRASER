annotate_strand_by_seq <- function(fds, genome="hg19", ...){
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
    stopifnot(is(genome, "BSgenome"))


    # get ranges of interest
    gr <- granges(granges(fds))

    # get the dinucleotide seq
    gr_start <- resize(gr, 2, fix = "start")
    gr_end   <- resize(gr, 2, fix = "end")

    mcols(fds, type="j")[["psi5_dinuc"]] <- getSeq(genome, gr_start, as.character=TRUE)
    mcols(fds, type="j")[["psi3_dinuc"]] <- getSeq(genome, gr_end,   as.character=TRUE)

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
    ans <- strand_table(fds, ...)

    strand_stat <- c(GTpos = sum(ans$GTpos), GCpos = sum(ans$GCpos), GTneg = sum(ans$GTneg), GCneg = sum(ans$GCneg),
                     half_pos = sum(ans$half_pos), half_neg = sum(ans$half_neg), posneg = sum(ans$posneg), none = sum(ans$none), total = length(fds)
    )
    return(strand_stat)
}

strand_table <- function(fds, ...){

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

    return(data.table(GTpos, GCpos, GTneg, GCneg, pos, neg, posneg, half_pos, half_neg, none))
}


plot_reads_strand <- function(fds, type, strand_type){
    dtstrand <- strand_table(fds)

    currentType(fds) <- type

    idx <- dtstrand[, strand_type, with =FALSE][[1]]
    dt_reads <- data.table(K = as.vector(K(fds)[idx, ] + 1), N = as.vector(N(fds)[idx, ] + 1))

    ggplot(dt_reads, aes(N, K)) + stat_binhex(aes(fill = log10(..count..))) +
        scale_x_log10() + scale_y_log10() +
        labs(title = paste0("Strand type \"", strand_type, "\""))
}


