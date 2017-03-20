#'
#' this is the example betabinom test
#'
#' testing and examples
#'
#' @author Julien Gagneur
#' @noRd
.juliens_betabinom_test <- function(){
    require(VGAM)

    ## simulate data
    n    = 100   ## number of samples
    size = 1000  ## number of trials
    prob = 0.2   ## prob of success per trial

    ## simulate coverage (denominator of psi)
    N = rbinom(n,size,prob)


    mu  = 0.9   ## mean psi
    rho = 0.2   ## some dispersion parameter, see ?rbetabinom

    ## simulate split read counts (numerator of psi)
    y  <- rbetabinom(n, size = N, prob = mu, rho = rho)


    ## set the first one as an outlier
    y[1] <- 1
    plot(N,y)

    ## fitting
    countMatrix <- cbind(y=y, o=N-y)
    fit <- vglm(cbind(y, N-y) ~ 1, betabinomial)
    co  <- Coef(fit)

    ## one-sided p-value (alternative = "less")
    pv  <- pbetabinom(y, N, co[["mu"]], co[["rho"]])
    plot(-log10(pv))
}

fraserNames <- data.table(
        readType = c("splitReads", "splitReads", "nonSplicedReads"),
        psiType = c("psi3", "psi5", "sitePSI"),
        pvalName = c("pvalue_psi3", "pvalue_psi5", "pvalue_sitePSI")
)

#'
#' this tests each splice type for all samples
#'
#' @noRd
.testPsiWithBetaBinomial <- function(dataset, internBPPARAM,
                    pvalFun=.betabinVglmTest){

    # check, that the object is stored as HDF5 array!
    if(!"DelayedArray" %in% is(assays(dataset@splitReads)[["rawCounts"]])){
        message(date(), ": The data is not stored in a HDF5Array. ",
            "To improve the performance we will store now ",
            "the data in HDF5 format.")
        saveFraseRDataSet(dataset)
    }

    # test all 3 different types
    for(idx in 1:nrow(fraserNames)){
        pvalName <- fraserNames[idx,pvalName]
        pvals <- .testPsiWithBetaBinomialPerType(dataset,
                fraserNames[idx,readType], fraserNames[idx,psiType], pvalFun
        )
        h5File <- .getFraseRHDF5File(dataset, fraserNames[idx,readType], TRUE)
        assays(slot(dataset, fraserNames[idx,readType]))[[pvalName]] <-
            .addAssayAsHDF5(pvals, pvalName, h5File)

        gc()
    }

    # return the new datasets
    return(dataset)
}

#'
#' extracts raw counts and saves it as a temp hdf5 assay
#'
#' @noRd
.getCountsAsHDF5 <- function(dataset, readType, countType){

    # get count data
    counts <- .getAssayAsDataTable(
            slot(dataset, readType), countType
    )

    # save it as hfd5
    assayname <- paste0(readType, "_", countType)
    tmpfile <- file.path(outputFolder(dataset@settings),
            "cache", "tmp", "counts.h5")
    if(!dir.exists(dirname(tmpfile))){
        dir.create(dirname(tmpfile), recursive=TRUE)
    } else if(file.exists(tmpfile) &&
                assayname %in% h5ls(tmpfile)$name){
        unlink(tmpfile)
    }
    h5obj <- writeHDF5Array(as.matrix(counts), tmpfile, assayname)

    return(h5obj)
}


#'
#' calculates the pvalue per type (psi3,psi5,spliceSite) with beta-binomial
#'
#' @noRd
.testPsiWithBetaBinomialPerType <- function(dataset, readType, psiType,
                    pvalFun){

    message(date(), ": Calculate P-values for the ",
            psiType, " splice type ..."
    )

    # go over each group but no NA's
    group         <- sampleGroup(dataset@settings)
    addNoise      <- TRUE
    removeHighLow <- 0
    #removeHighLow <- length(group)/4

    # reads to test for abberent splicing (eg: nonSplicedReads)
    rawCounts <- .getCountsAsHDF5(dataset, readType, countType="rawCounts")

    # other reads (eg: splitReads)
    rawOtherCounts <- .getCountsAsHDF5(dataset, readType,
            paste0("rawOtherCounts_", psiType)
    )

    # swap rawCounts with others for intron retention
    if(readType == "nonSplicedReads"){
        tmp <- rawCounts
        rawCounts <- rawOtherCounts
        rawOtherCounts <- tmp
    }

    # select junctions to test take only junctions
    # with at least 5 reads in at least one sample
    toTest <- which(apply(rawCounts,1,max) >= 5)
    message(date(), ": Sites to test: ", length(toTest))

    # TODO how to group groups?
    pvalues_ls <- bplapply(toTest, rawCounts=rawCounts,
            rawOtherCounts=rawOtherCounts,
            BPPARAM=dataset@settings@parallel,
            addNoise=addNoise, removeHighLow=removeHighLow,
            FUN= function(idx, rawCounts, rawOtherCounts,
                    addNoise, removeHighLow){
        suppressPackageStartupMessages(require(FraseR))

        ## simulate split read counts (numerator of psi)
        y <- as.integer(as.vector(unlist(rawCounts[idx,])))

        ## simulate coverage (denominator of psi)
        N <- y + as.integer(as.vector(unlist(rawOtherCounts[idx,])))

        # TODO betabinom fails if only zeros in one row
        # add noise to avoid zero variants in alternative reads
        if(addNoise){
            N <- N + sample(c(0,1),dim(rawCounts)[2], replace = TRUE)
        }

        # count matrix per site
        countMatrix <- cbind(y=y, o=N-y)

        # remove outliers (highest and lowest sample)
        tmp_removeHighLow = removeHighLow
        while(tmp_removeHighLow > 0){
            tmp_removeHighLow = tmp_removeHighLow -1
            countMatrix <- countMatrix[!1:dim(countMatrix)[1] %in%
                c(
                    which.min(rowSums(countMatrix)),
                    which.max(rowSums(countMatrix))
                ),
            ]
        }

        # plot(log10(N+1),log10(y+1))

        ## fitting
        pv_res <- lapply(list(countMatrix),
                FUN=.tryCatchFactory(pvalFun), y=y, N=N
        )[[1]]

        #
        # put pvalues into correct boundaries
        if(is.null(pv_res[[1]])){
            pv_res[[1]] <- rep(as.numeric(NA), dim(rawCounts)[2])
        }
        pv_res[[1]][which(pv_res[[1]] > 1)] <- as.numeric(NA)

        # plot(-log10(na.omit(pv_res[[1]])))
        return(pv_res)
    })

    # warnings
    warntable <- sort(table(gsub("^\\d+ diagonal ele", "xxx diagonal ele",
                                 unlist(sapply(pvalues_ls, "[[", "warn"))
    )))
    if(length(warntable) > 0){
        message(paste(collapse = "\n", c(date(),
            "Warnings in VGLM code while computing pvalues:\n",
            sapply(1:length(warntable), function(idx) paste(
                "\t", warntable[idx], "x", names(warntable)[idx]
            ))
        )))
    }

    # errors
    errotable <- sort(table(gsub("^\\d+ diagonal ele", "xxx diagonal ele",
                                 unlist(sapply(pvalues_ls, "[[", "err"))
    )))
    if(length(errotable) > 0){
        message(paste(collapse = "\n", c(
            "\nErrors in VGLM code while computing pvalues:\n",
            sapply(1:length(errotable), function(idx) paste(
                "\t", errotable[idx], "x", names(errotable)[idx]
            ))
        )))
    }

    # extract pvalues
    pvalues <- do.call(rbind,lapply(pvalues_ls, "[[", 1))
    pvalues_full <- matrix(as.numeric(NULL),
            nrow=dim(rawCounts)[1], ncol=dim(rawCounts)[2]
    )
    pvalues_full[toTest,] <- pvalues

    # transform it to a DataFrame and return it
    return(.asDataFrame(pvalues_full, dataset@settings@sampleData[,sampleID]))
}

#'
#' calculate the pvalues with vglm and the betabinomial functions
#'
#' @noRd
.betabinVglmTest <- function(countMatrix, y, N){
    fit <- vglm(countMatrix ~ 1, betabinomial)
    co  <- Coef(fit)

    ## one-sided p-value (alternative = "less")
    pbetabinom(y, N, co[["mu"]], co[["rho"]])
}


#'
#' calculate the pvalues with method of moments and the betabinomial functions
#' https://en.wikipedia.org/wiki/Beta-binomial_distribution#Method_of_moments
#'
#' @noRd
.betabinMMTest <- function(mat, y, N){
    # estimate mu1 and mu2
    Ns <- dim(mat)[1]
    n  <- ceiling(mean(rowSums(mat)))
    m1 <- 1/Ns * sum(mat[,"y"])
    m2 <- 1/Ns * sum(mat[,"y"]^2)

    # estimate the shapes
    nominator <- n*(m2/m1 - m1 - 1) + m1
    a <- (n*m1 - m2) / nominator
    b <- (n - m1)*(n - m2/m1) / nominator

    # return the probability
    pbetabinom.ab(y, N, max(0.1, a), max(0.1, b))
}

testing <- function(){
    mat <- countMatrix
    message("N:  ", Ns, "\nn:  ", n, "\nm1:  ", m1, "\nm2:  ", m2, "\na:  ", a, "\nb:  ", b)

    pF(10,n,a,b)
    plot(dbetabinom.ab(1:12, 12, a, b))
    cm <- matrix(ncol=2,
            c(0:12, c(3, 24, 104, 286, 670, 1033, 1343, 1112, 829, 478, 181, 45, 7))
    )
    colnames(cm) <- c("x", "y")

    fit <- vglm(countMatrix ~ 1, betabinomial)
    co  <- Coef(fit)

    ## one-sided p-value (alternative = "less")
    pbetabinom(y, N, co[["mu"]], co[["rho"]])

}

##
## error/warning catching functions by martin morgan
## http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
##
.tryCatchFactory <- function(fun){
    function(...) {
        warn <- err <- NULL
        res <- withCallingHandlers(
            tryCatch(fun(...), error=function(e) {
                err <<- conditionMessage(e)
                NULL
            }), warning=function(w) {
                warn <<- append(warn, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        list(res, warn=warn, err=err)
    }
}
