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
        dataset <- saveFraseRDataSet(dataset)
    }

    # test all 3 different types
    for(idx in 1:nrow(fraserNames)){
        pvalName <- fraserNames[idx,pvalName]
        readType <- fraserNames[idx,readType]
        psiType  <- fraserNames[idx,psiType]
        dataset <- .testPsiWithBetaBinomialPerType(
                dataset, readType, psiType, pvalName, pvalFun
        )
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
                    pvalname, pvalFun){

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
        # idx <- toTest[[2]]
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
            pv_res[[1]] <- list(
                pval = rep(as.numeric(NA), dim(rawCounts)[2]),
                result = paste0("ERROR: ", pv_res$err)
            )
        }
        pv_res[[1]]$pval[which(pv_res[[1]]$pval > 1)] <- as.numeric(NA)

        # plot(-log10(na.omit(pv_res[[1]])))
        return(pv_res)
    })

    # print the vglm infos
    for(type in c("warn", "err")){
        infoChar <- .vglmInfos2character(pvalues_ls, type)
        if(infoChar != ""){
            message(infoChar)
        }
    }

    # extract additional informations
    mcol_ls <- list(
        alpha   = as.vector(sapply(pvalues_ls, function(x) x[[1]]$alpha)),
        beta    = as.vector(sapply(pvalues_ls, function(x) x[[1]]$beta)),
        timing  = as.vector(sapply(pvalues_ls, function(x) x[[1]]$timing)),
        errStr  = .vglmInfos2character(pvalues_ls, "err"),
        warnStr = .vglmInfos2character(pvalues_ls, "warn")
    )
    for(i in seq_along(mcol_ls)){
        name <- paste0(psiType, "_", names(mcol_ls)[[i]])
        mcols(slot(dataset, readType))[[name]] <- NA
        mcols(slot(dataset, readType))[[name]][toTest] <- mcol_ls[[i]]
    }

    # extract pvalues
    pvalues <- do.call(rbind,lapply(pvalues_ls, function(x) x[[1]]$pval))
    pvalues_full <- matrix(as.numeric(NULL),
            nrow=dim(rawCounts)[1], ncol=dim(rawCounts)[2]
    )
    pvalues_full[toTest,] <- pvalues

    # transform it to a hdf5 assay and save it to the dataset
    h5File  <- .getFraseRHDF5File(dataset, readType, TRUE)
    pval_df <- .asDataFrame(pvalues_full, samples(dataset@settings))
    h5Assay <- .addAssayAsHDF5(pval_df, pvalname, h5File)
    assays(slot(dataset, readType))[[pvalname]] <- h5Assay

    # save the pvalue calculations in the SE ranged object
    return(dataset)
}

#'
#' calculate the pvalues with vglm and the betabinomial functions
#'
#' @noRd
.betabinVglmTest <- function(countMatrix, y, N){
    time <- system.time({
        # get fit
        fit <- vglm(countMatrix ~ 1, betabinomial)

        # get the shape values
        co  <- Coef(fit)
        prob <- co[["mu"]]
        rho  <- co[["rho"]]
        alpha <- prob * (1 - rho)/rho
        beta  <- (1 - prob) * (1 - rho)/rho

        # get the pvalues
        # one-sided p-value (alternative = "less")
        pval <- pbetabinom.ab(y, N, alpha, beta)
    })["elapsed"]

    return(list(
        pval = pval,
        alpha = alpha,
        beta = beta,
        timing = time
    ))
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

#'
#' error/warning catching functions by martin morgan
#' http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#'
#' @noRd
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


#'
#' formate the VGLM error/warning output in a nice readable way
#'
#' @noRd
.vglmInfos2character <- function(res, type){
    infoList <- unlist(sapply(res, "[[", type))
    infoList <- gsub("^\\d+ diagonal ele", "xxx diagonal ele", infoList)
    infoTable <- sort(table(infoList))

    if(length(infoTable) == 0){
        return("")
    }

    return(paste(collapse = "\n", c(paste0(date(), ": ",
            ifelse(type == "warn", "Warnings", "Errors"),
            " in VGLM code while computing pvalues:"),
            .table2character(warntable)
    )))
}


#'
#' convert a table to a string like the output of table
#'
#' @noRd
.table2character <- function(table){
    sapply(1:length(table), function(idx) paste(
        "\t", table[idx], "x", names(table)[idx]
    ))
}



.testingbetabinom <- function(){
    fds <- createTestFraseRDataSet()

    source("FraseR/R/beta-binomial-testing.R")
    source("FraseR/R/helper-functions.R")
    pvalFun=.betabinVglmTest
    idx <- 1
    readType <- fraserNames[idx,readType]
    psiType <- fraserNames[idx,psiType]
    internBPPARAM <- MulticoreParam(4,progressbar=TRUE)
    parallel(fds@settings) <- MulticoreParam(4,progressbar=TRUE)
    dataset <- fds

    ## one-sided p-value (alternative = "less")
    pbetabinom(y, N, co[["mu"]], co[["rho"]])

}

