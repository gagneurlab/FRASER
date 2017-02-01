
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
y = rbetabinom(n, size = N, prob = mu, rho = rho) 


## set the first one as an outlier
y[1] = 1 
plot(N,y)

## fitting
fit <- vglm(cbind(y, N-y) ~ 1, betabinomial)
co = Coef(fit)
## one-sided p-value (alternative = "less") 
pv = pbetabinom(y, N, co[["mu"]], co[["rho"]])
plot(-log10(pv))

readType <- "splitReads"
psiType  <- "psi3"

index = c(1:dim(rawCounts)[1])
timmdc1_index <- c(45,48,52,56,58,60)
mcoln1_index  <- c(102,104,107,109)
#  index <- timmdc1_index
rawCounts[index,]
# idx = 1
#    c1 = rawCounts[index,]
#   c2 = rawOtherCounts[index,]

params <- dataset@settings@parallel
bpcatchErrors(params) <- TRUE
bpstopOnError(params) <- FALSE
params <- SerialParam()

#'
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
    group   <- sampleGroup(dataset@settings)
    
    # reads to test for abberent splicing (eg: nonSplicedReads)
    rawCounts <- FraseR:::.getAssayAsDataTable(slot(dataset, readType), "rawCounts")
    
    # other reads (eg: splitReads)
    rawOtherCounts <- FraseR:::.getAssayAsDataTable(slot(dataset, readType), paste0("rawOtherCounts_", psiType))
  
    # select junctions to test
    toTest <- which(apply(rawCounts,1,max) >= 5)
    
    # TODO how to group groups?
    pvalues_ls <- bplapply(toTest, rawCounts=rawCounts, 
            rawOtherCounts=rawOtherCounts,
            BPPARAM=params,
            FUN= function(idx, rawCounts, rawOtherCounts){
        library(FraseR)
        require(VGAM)
        
        ## simulate split read counts (numerator of psi)
        y <- as.integer(as.vector(unlist(rawCounts[idx,])))
       
        ## simulate coverage (denominator of psi)
        N <- (y + as.integer(unlist(rawOtherCounts[idx,]))
             # add noise to avoid zero variants
             + sample(c(0,1),dim(rawCounts)[2], replace = T)
        )
        
        # count matrix per site
        countMatrix <- cbind(y=y, o=N-y)
        
        # remove outliers (highest and lowest sample)
        rmMaxNMin <- 1:dim(countMatrix)[1] %in% c(
            which.min(rowSums(countMatrix)), 
            which.max(rowSums(countMatrix))
        )
        
        
        # print(cbind(countMatrix, rmMaxNMin))
        # plot(N,y)

        ## fitting
        pv <- rep(as.numeric(NA), dim(rawCounts)[2])
        tryCatch({
            fit <- vglm(countMatrix[!rmMaxNMin,] ~ 1, betabinomial)
            co  <- Coef(fit)
            
            ## one-sided p-value (alternative = "less")
            pv  <- pbetabinom(y, N, co[["mu"]], co[["rho"]])
        }, warning = function(w) str(w), error = function(e) str(e))
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


#'
#' calculates the pvalue per group with fisher
#' 
#' @noRd
.testPsiWithFisherPerGroup <- function(dataset, groupID, rawCounts, rawOtherCounts, internBPPARAM){
    # get group to test
    group <- sampleGroup(dataset@settings)
    group2Test <- group == groupID
    group2Test[is.na(group2Test)] <- FALSE
    
    fullFisherTable <- data.table(
        TP=rowSums(rawCounts[     , group2Test,with=FALSE]),
        FP=rowSums(rawCounts[     ,!group2Test,with=FALSE]),
        FN=rowSums(rawOtherCounts[, group2Test,with=FALSE]),
        TN=rowSums(rawOtherCounts[,!group2Test,with=FALSE])
    )
    
    # test only where at least the group has one read
    fisherTableToTest <- fullFisherTable[TP+FN > 0]
    pvalues <- unlist(bplapply(1:nrow(fisherTableToTest), BPPARAM = internBPPARAM, fisherTableToTest=fisherTableToTest,
                               function(idx, fisherTableToTest){
                                   fisher.test(matrix(as.integer(fisherTableToTest[idx]), nrow=2))$p.value 
                               }
    ))
    
    # add NAs wher the test group did not had any read
    fullFisherTable[,pvalue:=as.numeric(NA)]
    fullFisherTable[TP+FN>0,pvalue:=pvalues]
    return(fullFisherTable[,pvalue])
}





