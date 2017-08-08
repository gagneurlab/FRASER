#'
#' calculates the pvalue per type (psi3,psi5,spliceSite) with beta-binomial
#'
#' @noRd
pvalueByBetaBinomialPerType <- function(fds, aname, psiType, pvalFun,
            minCov=5, alternative="less"){

    message(date(), ": Calculate P-values for the ",
            psiType, " splice type ...")

    # go over each group but no NA's
    group         <- condition(fds)
    addNoise      <- TRUE

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
    pvalues_ls <- bplapply(toTest, rawCounts=rawCounts,
            rawOtherCounts=rawOtherCounts, pvalFun,
            BPPARAM=parallel(fds),
            addNoise=addNoise, removeHighLow=removeHighLow,
            FUN= function(idx, rawCounts, rawOtherCounts,
                    addNoise, removeHighLow, pvalFun){
        # idx <- toTest[[2]]
        suppressPackageStartupMessages(require(FraseR))

        ## get read coverage (y == reads of interests, N == all reads)
        y <- as.integer(as.matrix(rawCounts[idx,]))
        N <- y + as.integer(as.matrix(rawOtherCounts[idx,]))

        # TODO betabinom fails if only zeros in one row
        # add noise to avoid zero variants in alternative reads
        if(sd(N) == 0 & addNoise){
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
        timing <- sum(system.time(gcFirst=FALSE, expr={
            pv_res <- lapply(list(countMatrix), FUN=tryCatchFactory(pvalFun),
                    y=y, N=N, alternative=alternative)[[1]]
        })[c("user.self", "sys.self")])

        #
        # put pvalues into correct boundaries
        if(is.null(pv_res[[1]])){
            pv_res[[1]] <- list(
                pval = rep(as.numeric(NA), dim(rawCounts)[2]),
                alpha = NA,
                beta = NA
            )
        }
        pv_res[[1]]$pval[which(pv_res[[1]]$pval > 1)] <- as.numeric(NA)

        # add timing
        pv_res[[1]]$timing <- timing

        # add index for later sorting
        pv_res[[1]]$idx <- idx

        return(pv_res)
    })

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
betabinVglmTest <- function(cMat, alternativ="less", y=cMat[,1], N=cMat[,1] + cMat[,2]){
    # get fit
    fit <- vglm(cMat ~ 1, betabinomial)

    # get the shape values
    co  <- Coef(fit)
    prob <- co[["mu"]]
    rho  <- co[["rho"]]
    alpha <- prob * (1 - rho)/rho
    beta  <- (1 - prob) * (1 - rho)/rho

    # get the pvalues
    # one-sided p-value (alternative = "less")
    pval <- pbetabinom.ab(y, N, alpha, beta)

    # two sieded test
    if(alternative == "two.sided"){
        pval <- apply(cbind(pval, 1-pval), 1, min)*2
    # one-sided greater test
    } else if(alternative == "greater"){
        pval <- 1-pval
    }

    return(list(
        pval = pval,
        alpha = alpha,
        beta = beta
    ))
}


#'
#' error/warning catching functions by martin morgan
#' http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
#'
#' @noRd
tryCatchFactory <- function(fun){
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


