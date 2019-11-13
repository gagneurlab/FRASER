
#'
#' calculates the pvalue per type (psi3,psi5,spliceSite) with beta-binomial
#'
#' @examples
#' fds <- countRNAData(createTestFraseRSettings())
#' fds <- calculatePSIValues(fds)
#' debug(FraseR:::betabinVglmTest)
#' fds <- pvalueByBetaBinomialPerType(fds, "testPval", "psi5",
#'         FraseR:::betabinVglmTest)
#'
#' @noRd
pvalueByBetaBinomialPerType <- function(fds, aname, psiType, pvalFun,
            minCov=5, alternative="two.sided", timeout=300,
            BPPARAM=bpparam(), returnFit=FALSE){

    message(date(), ": Calculate P-values for the ",
            psiType, " splice type ...")

    # go over each group but no NA's
    group         <- condition(fds)
    addNoise      <- FALSE

    # raw reads to test for abberent splicing
    rawCounts      <- counts(fds, type=psiType, side="ofInterest")
    rawOtherCounts <- counts(fds, type=psiType, side="otherSide")

    # swap rawCounts with others for intron retention
    if(psiType == "psiSite"){
        tmp <- rawCounts
        rawCounts <- rawOtherCounts
        rawOtherCounts <- tmp
        rm(tmp)
    }

    # select junctions to test take only junctions
    # with at least 5 reads in at least one sample
    # toTest <- which(apply(rawCounts,1,median) >= 5)
    toTest <- which(apply(rawCounts,1,max) >= minCov)
    message(date(), ": Sites to test: ", length(toTest))

    # TODO: IMPORTANT:
    # Tests are not equal in there runtime and
    # long running test tent to cluster at genomic positions
    # therefore shuffel positions to remove this runtime bias
    toTest <- sample(toTest)

    # TODO how to group groups?
    testFUN <- function(idx, rawCounts, rawOtherCounts, addNoise, pvalFun,
                        timeout, returnFit=FALSE, ...){

        ## get read coverage (y == reads of interests, N == all reads)
        y <- as.integer(as.matrix(rawCounts[idx,]))
        N <- y + as.integer(as.matrix(rawOtherCounts[idx,]))

        # TODO betabinom fails if only zeros in one row
        # add noise to avoid zero variants in alternative reads
        if(sd(N) == 0 | sd(N-y) == 0 | addNoise == TRUE){
            N <- N + sample(c(0,1), ncol(rawCounts), replace=TRUE)
        }

        # count matrix per site
        # plot(log10(N+1),log10(y+1))
        countMatrix <- cbind(y=y, o=N-y)

        ## fitting
        timing <- sum(system.time(gcFirst=FALSE, expr={
            pv_res <- lapply(list(countMatrix), y=y, N=N, ...,
                    FUN=tryCatchFactory(pvalFun, timeout=timeout))[[1]]
        })[c("user.self", "sys.self")])

        #
        # put pvalues into correct boundaries
        if(is.null(pv_res[[1]])){
            pv_res[[1]] <- list(
                pval = rep(as.numeric(NA), ncol(rawCounts)),
                alpha = NA,
                beta = NA
            )
        }
        pv_res[[1]]$pval[which(pv_res[[1]]$pval > 1)] <- as.numeric(NA)

        # add timing
        pv_res[[1]]$timing <- timing

        # add index for later sorting
        pv_res[[1]]$idx <- idx

        if(isFALSE(returnFit)){
            pv_res[[1]]$fit <- NULL
        }
        return(pv_res)
    }

    FUN <- function(...){ testFUN(...) }
    gc()

    pvalues_ls <- bplapply(toTest, rawCounts=rawCounts, alternative=alternative,
                    rawOtherCounts=rawOtherCounts, pvalFun, timeout=timeout,
                    returnFit=returnFit, addNoise=addNoise, FUN=FUN,
                    BPPARAM=BPPARAM)

    # sort table again
    pvalues_ls <- pvalues_ls[order(toTest)]
    toTest <- sort(toTest)

    # print the vglm infos
    for(type in c("warn", "err")){
        infoChar <- vglmInfos2character(pvalues_ls, type)
        if(infoChar != ""){
            message(infoChar)
        }
    }

    # extract additional informations
    mcol_ls <- list(
        tested  = TRUE,
        alpha   = as.vector(sapply(pvalues_ls, function(x) x[[1]]$alpha)),
        beta    = as.vector(sapply(pvalues_ls, function(x) x[[1]]$beta)),
        timing  = as.vector(sapply(pvalues_ls, function(x) x[[1]]$timing)),
        idx     = as.vector(sapply(pvalues_ls, function(x) x[[1]]$idx)),
        errStr  = vglmInfos2character(pvalues_ls, "err"),
        warnStr = vglmInfos2character(pvalues_ls, "warn")
    )
    for(i in seq_along(mcol_ls)){
        name <- paste0(psiType, "_", names(mcol_ls)[[i]])
        mcols(fds, type=psiType)[[name]] <- as(NA, class(mcol_ls[[i]]))
        mcols(fds, type=psiType)[[name]][toTest] <- mcol_ls[[i]]
    }

    # extract pvalues
    pvalues <- do.call(rbind,lapply(pvalues_ls, function(x) x[[1]]$pval))
    pvalues_full <- matrix(as.numeric(NA),
            nrow=dim(rawCounts)[1], ncol=dim(rawCounts)[2]
    )
    pvalues_full[toTest,] <- pvalues

    # transform it to a hdf5 assay and save it to the dataset
    assays(fds, type=psiType)[[aname]] <- pvalues_full

    # save the pvalue calculations in the SE ranged object
    return(fds)
}

#'
#' calculate the pvalues with vglm and the betabinomial functions
#'
#' @noRd
betabinVglmTest <- function(cMat, alternative="two.sided",
                    y=cMat[,1], N=cMat[,1] + cMat[,2]){
    # get fit
    fit <- vglm(cMat ~ 1, betabinomial)

    # get the shape values
    co  <- Coef(fit)
    prob <- co[["mu"]]
    rho  <- co[["rho"]]
    alpha <- prob * (1 - rho)/rho
    beta  <- (1 - prob) * (1 - rho)/rho

    # get the pvalues only for non zero values
    naValues <- is.na(y + N) | y + N == 0

    # one-sided p-value (alternative = "less")
    pval <- pbetabinom.ab(y[!naValues], N[!naValues], alpha, beta)
    dval <- dbetabinom.ab(y[!naValues], N[!naValues], alpha, beta)
    pval <- pmin(pval, 1)
    dval <- pmin(dval, 1)

    # two sieded test
    if(startsWith("two.sided", alternative)){
        pval <- pmin(pval * 2, (1 - pval + dval) * 2, 1)
    # one-sided greater test
    } else if(startsWith("greater", alternative)){
        pval <- pmin(1 - pval + dval, 1)
    } else if(!startsWith("less", alternative)){
        stop("")
    }

    # set NA values for non tested samples
    if(any(naValues)){
        tmp <- rep(NA, length(naValues))
        tmp[!naValues] <- pval
        pval <- tmp
    }

    return(list(
        pval  = pval,
        alpha = alpha,
        beta  = beta,
        prob  = prob,
        rho   = rho,
        fit   = fit
    ))
}


#'
#' error/warning catching functions by martin morgan
#' http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#'
#' @noRd
tryCatchFactory <- function(fun, timeout=300){
    function(...) {
        warn <- err <- NULL
        res <- withCallingHandlers(
            tryCatch(withTimeout(fun(...), timeout=timeout),
                error=function(e) {
                    err <<- conditionMessage(e)
                    NULL
                }),
            warning=function(w) {
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
vglmInfos2character <- function(res, type){
    infoList <- unlist(sapply(res, "[[", type))
    infoList <- gsub("^\\d+ diagonal ele", "xxx diagonal ele", infoList)
    infoTable <- sort(table(infoList))

    if(length(infoTable) == 0){
        return("")
    }

    return(paste(collapse = "\n", c(paste0(date(), ": ",
            ifelse(type == "warn", "Warnings", "Errors"),
            " in VGLM code while computing pvalues:"),
            table2character(infoTable)
    )))
}


#'
#' convert a table to a string like the output of table
#'
#' @noRd
table2character <- function(table){
    sapply(1:length(table), function(idx) paste(
        "\t", table[idx], "x", names(table)[idx]
    ))
}


testing <- function(){
    ocMat <- cMat

    ccMat <- ocMat
    rcMat <- rowSums(ccMat)
    targetQ <- 0.75
    targetDepth <- quantile(rowSums(ccMat), targetQ)

    sizeFactor <- targetDepth / rcMat

    ncMat <- ceiling(ccMat * sizeFactor)+1

    cMat <- ncMat[2:101,]
    y <- cMat[,1]
    N <- rowSums(cMat)

    plot(-log10((1:length(pval))/length(pval)), -log10(sort(pval)))
    abline(0,1)
    grid()
}


#' @export
cpp_truncNll_db <- function(par, H, k, n, rho){
    truncNLL_db(par=par, H=H, k=k, n=n, rho=rho)
}

#' @export
cpp_truncGrad_db <- function(par, H, k, n, rho){
    truncGrad_db(par=par, H=H, k=k, n=n, rho=rho)
}
