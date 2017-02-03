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
    fit <- vglm(cbind(y, N-y) ~ 1, betabinomial)
    co  <- Coef(fit)
    
    ## one-sided p-value (alternative = "less") 
    pv  <- pbetabinom(y, N, co[["mu"]], co[["rho"]])
    plot(-log10(pv))
}

#'
#' this tests each splice type for all samples
#' 
#' @noRd
.testPsiWithBetaBinomial <- function(dataset){
    
    # test all 3 different types
    assays(dataset@splitReads)$pvalue_psi3 <- 
        .testPsiWithBetaBinomialPerType(dataset, "splitReads", "psi3")
    assays(dataset@splitReads)$pvalue_psi5 <- 
        .testPsiWithBetaBinomialPerType(dataset, "splitReads", "psi5")
    assays(dataset@nonSplicedReads)$pvalue_sitePSI <- 
        .testPsiWithBetaBinomialPerType(dataset, "nonSplicedReads", "sitePSI")
    
    # return the new datasets
    return(dataset)
}

#'
#' calculates the pvalue per type (psi3,psi5,spliceSite) with beta-binomial
#' 
#' @noRd
.testPsiWithBetaBinomialPerType <- function(dataset, readType, psiType){
    # go over each group but no NA's
    group         <- sampleGroup(dataset@settings)
    addNoise      <- TRUE
    removeHighLow <- length(group)/4
    
    # reads to test for abberent splicing (eg: nonSplicedReads)
    rawCounts <- FraseR:::.getAssayAsDataTable(slot(dataset, readType), "rawCounts")
    
    # other reads (eg: splitReads)
    rawOtherCounts <- FraseR:::.getAssayAsDataTable(slot(dataset, readType), paste0("rawOtherCounts_", psiType))
  
    # swap rawCounts with others for intron retention
    if(readType == "nonSplicedReads"){
        tmp <- rawCounts
        rawCounts <- rawOtherCounts
        rawOtherCounts <- tmp
    }
    
    # select junctions to test take only junctions 
    # with at least 5 reads in at least one sample
    toTest <- which(apply(rawCounts,1,max) >= 5)
    
    # TODO how to group groups?
    pvalues_ls <- bplapply(toTest, rawCounts=rawCounts, 
            rawOtherCounts=rawOtherCounts,
            BPPARAM=dataset@settings@parallel,
            addNoise=addNoise, removeHighLow=removeHighLow,
            FUN= function(idx, rawCounts, rawOtherCounts, 
                          addNoise, removeHighLow){
        library(FraseR)
        require(VGAM)
        
        ## simulate split read counts (numerator of psi)
        y <- as.integer(as.vector(unlist(rawCounts[idx,])))
       
        ## simulate coverage (denominator of psi)
        N <- y + as.integer(unlist(rawOtherCounts[idx,]))
        
        # TODO betabinom fails if only zeros in one row
        # add noise to avoid zero variants in alternative reads
        if(addNoise){
            N <- N + sample(c(0,1),dim(rawCounts)[2], replace = T)
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
        
        # bootstrap sample data if sample count to low
        if(dim(countMatrix)[1] < 25){
            warning(paste0("Sample count (", 
                    dim(countMatrix)[1], 
                    ") to low. Pseudo bootstrap is used to have 25 samples."
            ))   
            samples_to_take <- c(
                1:dim(countMatrix)[1], 
                sample(1:dim(countMatrix)[1], 
                       max(0, 25 - dim(countMatrix)[1]), 
                       replace = TRUE
                )
            )
            countMatrix <- countMatrix[samples_to_take,]
        }
        
        # plot(log10(N+1),log10(y+1))

        ## fitting
        pv <- rep(as.numeric(NA), dim(rawCounts)[2])
        tryCatch({
            fit <- vglm(countMatrix ~ 1, betabinomial)
            co  <- Coef(fit)
            
            ## one-sided p-value (alternative = "less")
            pv  <- pbetabinom(y, N, co[["mu"]], co[["rho"]])
        }, error = function(e) str(e))
        pv[which(pv > 1)] <- as.numeric(NA)
        
        # plot(-log10(na.omit(pv)))
        return(pv)
    })
    
    
    pvalues <- do.call(rbind,pvalues_ls)
    
    pvalues_full <- matrix(as.numeric(NULL), nrow=dim(rawCounts)[1], ncol=dim(rawCounts)[2])
    pvalues_full[toTest,] <- pvalues
    
    # transform it to a DataFrame and return it
    return(FraseR:::.asDataFrame(pvalues_full, dataset@settings@sampleData[,sampleID]))
}
