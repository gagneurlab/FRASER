context("Gradient testing")

test_that("CPP gradient", {
    testFun <- function(){
        getX <- function(K, N){
            t(qlogis((K + 0.5)/(N + 1)))
        }
        
        numericLossGrad <- function(fn, epsilon, w, ..., BPPARAM=bpparam()){
            grad <- unlist(bplapply(seq_along(w), fn=fn, w=w, epsilon=epsilon,
                            ..., BPPARAM=BPPARAM, function(i, fn, epsilon, w, 
                                            ...){
                eps <- integer(length(w))
                eps[i] <- epsilon
                (fn(w + eps, ...) - fn(w - eps, ...))/(2*epsilon)
            }))
            return(grad)
        }
        
        
        #
        # Simulations
        #
        # some expample matrices we, wd, b, rho, K, N and x to test the gradient
        n = 100
        j = 300
        q = 3
        
        we   <- matrix(runif(j*q, -5, 5), nrow=j)
        wd   <- matrix(runif(j*q, -5, 5), nrow=j)
        b    <- runif(j, -1, 1)
        rho  <- matrix(runif(j, 0, 1), ncol=n, nrow=j)
        prob <- matrix(runif(j, 0, 1), ncol=n, nrow=j)
        
        cov <- matrix(ceiling(rlnorm(j, 5)) + 10, ncol=n, nrow=j)
        disp <- rlnorm(j, 4) + 1
        
        N <- matrix(rnbinom(n*j, mu=cov, size=rep(disp, each=n)), 
                ncol=n, byrow=FALSE)
        K <- matrix(rbetabinom(n*j, size=N, prob=prob, rho=rho),
                ncol=n, byrow=FALSE)
        
        #
        # Gradient and loss for D and b (is calculated per gene)
        #
        mybprange <- seq_len(j)
        res <- bplapply(mybprange, K=K, N=N, wd=wd, we=we, b=b, rho=rho, 
                        function(idx, K, N, rho, b, wd, we){
            message(idx)
            
            maxY <- 700
            minU <- 0 #1e-20
        
            ki   <- K[idx,]
            ni   <- N[idx,]
            x    <- getX(K, N)
            H    <- x %*% we
            rhoi <- rho[idx]
            bi   <- b[idx]
            wdi  <- wd[idx,]
            par <- c(bi, wdi)
        
            estimateLgamma_alpha <- function(x, r){
        
                fixed <- lgamma(plogis(-30) * r)
                # step <- lgamma(plogis(-30) * -r) - lgamma(plogis(-29) * -r)
                est <- fixed - 1*(x+30)
        
                return(est)
            }
        
            estimateLgamma_beta <- function(x, r){
        
                fixed <- lgamma((plogis(30)-1) * -r)
                # step <- lgamma((plogis(30)-1) * -r) - 
                #         lgamma((plogis(29)-1) * -r)
                est <- fixed + 1*(x-30)
        
                return(est)
            }
        
            loosDb <- function(par, ki, ni, H, rhoi){
                wdi <- par[2:length(par)]
                bi <- par[1]
        
                yi <- H %*% wdi + bi
                # yi <- pmax(pmin(yi, maxY), -maxY)
                ey <- exp(yi)
        
                ar <- (1 - rhoi)/rhoi
                br <- (rhoi - 1)/rhoi
                ph <- ey/(1+ey)
                u  <- -1/(ey + 1)
        
                ph[!is.finite(ph)] <- vapply(yi[!is.finite(ph)],
                                            function(y){
                                                if(y < 0){0}else{1} },
                                            double(1))
                u[!is.finite(u)] <- vapply(yi[!is.finite(u)],
                                            function(y){
                                                if(y < 0){-1}else{0} },
                                            double(1))
        
                alpha  <- lgamma(ph * ar)
                alphaK <- lgamma(ph * ar + ki + 0.5)
                beta   <- lgamma(u * br)
                betaNK <- lgamma(u * br + ni - ki + 0.5)
        
                # est <- abs(yi)
                # alpha[is.infinite(alpha)] <- est[is.infinite(alpha)]
                # beta[is.infinite(beta)] <- est[is.infinite(beta)]
                alpha[is.infinite(alpha)] <- 
                    estimateLgamma_alpha(yi,ar)[is.infinite(alpha)]
                beta[is.infinite(beta)] <- 
                    estimateLgamma_beta(yi,ar)[is.infinite(beta)]
        
                mean(alpha + beta - alphaK - betaNK)
            }
        
            gradDb <- function(par, ki, ni, H, rhoi){
                wdi <- par[2:length(par)]
                bi  <- par[1]
        
                yi <- H %*% wdi + bi
                # yi <- pmax(pmin(yi, maxY), -maxY)
                ey <- exp(yi)
        
                ar <- (1 - rhoi)/rhoi
                br <- (rhoi - 1)/rhoi
                ph <- ey/(1+ey)
                u  <- -1/(ey + 1)
                v  <- ey/(1+ey)^2
        
                ph[!is.finite(ph)] <- vapply(yi[!is.finite(ph)],
                                            function(y){
                                                if(y < 0){0} else{1} },
                                            double(1))
                u[!is.finite(u)] <- vapply(yi[!is.finite(u)],
                                            function(y){
                                                if(y < 0){-1}else{0} },
                                            double(1))
        
        
        
                alpha  <- digamma(ph * ar) * ar * v
                alphaK <- digamma(ph * ar + ki + 0.5) * ar * v
                beta   <- digamma(u * br) * br * v
                betaNK <- digamma(u * br + ni - ki + 0.5) * br * v
        
                v[is.nan(v)] <- 0
                alpha[v == 0] <- vapply(yi[v == 0],
                                                function(y){
                                                    if(y < 0){-1}else{0} },
                                                double(1))
                alpha[!is.finite(alpha)] <- vapply(yi[!is.finite(alpha)],
                                                function(y){
                                                    if(y < 0){-1}else{0} },
                                                double(1))
                beta[v == 0] <- vapply(yi[v == 0],
                                        function(y){
                                            if(y < 0){0}else{1} },
                                        double(1))
                beta[!is.finite(beta)] <- vapply(yi[!is.finite(beta)],
                                                function(y){
                                                    if(y < 0){0}else{1} },
                                                double(1))
                alphaK[!is.finite(alphaK)] <- 0
                betaNK[!is.finite(betaNK)] <- 0
        
                idx2replace = ph == 0 | ph == 1
                # H[idx2replace,] <- 0
        
                grad_b <- mean(alpha + beta - alphaK - betaNK)
                grad_d_mat <- matrix((alpha + beta - alphaK - betaNK), nrow(H), 
                                    ncol(H)) * H
                #grad_d_mat[ph == 0 | ph == 1] <- vapply(yi[ph == 0 | ph == 1],
                #                                function(y){
                #                                  if(y < 0) {1} else { -1 } },
                #                                double(1))
        
                grad_d <- colMeans(grad_d_mat)
        
                grad_d_replace <- vapply(seq_along(grad_d), function(j){
                    badIdx = which(ph[j] == 0 | ph[j] == 1)
                    if(sum(yi[badIdx] < 0) < length(badIdx)){
                        -1
                    } else {
                        +1
                    }
                }, FUN.VALUE=numeric(1))
        
        
                #grad_d[idx2replace] <- grad_d_replace[idx2replace]
        
                c(grad_b, grad_d)
            }
        
            nlgv <- numericLossGrad(loosDb, 1e-8, par, ki=ki, ni=ni, H=H,
                    rhoi=rhoi)
            grv <- gradDb(par, ki, ni, H, rhoi)
            cppnlgv <- numericLossGrad(truncNLL_db, 1e-8, par, k=ki, n=ni, H=H, 
                    rho=rhoi, lambda = 0)
            cppgrv  <- as.vector(truncGrad_db(par=par, H=H, k=ki, 
                    n=ni, rho=rhoi, lambda = 0))
        
            list(nlgv=nlgv, grv=grv, cppnlgv=cppnlgv, cppgrv=cppgrv)
        })
        
        nlgv <- vapply(res, "[[", "nlgv", FUN.VALUE=numeric(ncol(wd)+1))
        grv  <- vapply(res, "[[", "grv", FUN.VALUE=numeric(ncol(wd)+1))
        cppnlgv <- vapply(res, "[[", "cppnlgv", FUN.VALUE=numeric(ncol(wd)+1))
        cppgrv  <- vapply(res, "[[", "cppgrv", FUN.VALUE=numeric(ncol(wd)+1))
        
        nlgv[,idx]
        grv[,idx]
        
        par(mfrow=c(1,2))
        plot(log10(sort(abs(nlgv - grv))))#, ylim=c(-11,1))
        plot(nlgv, grv); grid(); abline(0,1, col="red")
        plot(log10(sort(abs(cppnlgv - cppgrv)))) #, ylim=c(-11,1))
        plot(cppnlgv, cppgrv); grid(); abline(0,1, col="red")
        
        which(abs(nlgv - grv) > 0.01, arr.ind = TRUE)
        which(abs(nlgv - grv) > 10, arr.ind = TRUE)
        table(sign(nlgv) == sign(grv))
        which(!(sign(nlgv) == sign(grv)), arr.ind = TRUE)
        par(mfrow=c(1,1))
        
        plotInaccuracy <- function(idx, col, range=1e-7, xlegend=NULL,
                            ylegend=NULL){
        
            ki   <- K[idx,]
            ni   <- N[idx,]
            x    <- getX(K, N)
            H    <- x %*% we
            rhoi <- rho[idx]
            bi   <- b[idx]
            wdi  <- wd[idx,]
            par <- c(bi, wdi)
        
            if(is.null(ylegend)){
                l <- rep(0, q+1)
                l[col] <- range
                ylegend <- loosDb(par + l, ki, ni, H, rhoi)
            }
            if(is.null(xlegend)){
                xlegend <- -range
            }
        
            plot(seq(-range, range, length.out = 1000), 
                vapply(seq(-range, range, length.out = 1000),
                        function(x){
                            eps <- rep(0, q+1); eps[col] <- x; 
                            loosDb(par + eps, ki, ni, H, rhoi)}, 
                        double(1)),
                type = "l", xlab = "epsilon", ylab = "loss(D[q] + epsilon)", 
                main=paste("idx:", idx, "q:", col)); grid();
            abline(v=c(-1e-8, 1e-8), col = c("grey", "grey"), lty=c(2,2));
            abline(loosDb(par, ki, ni, H, rhoi), nlgv[col,idx], 
                    col="green", lty=2); 
            abline(loosDb(par, ki, ni, H, rhoi), grv[col,idx], 
                    col="red", lty=2);
            legend(xlegend, ylegend, legend=c("loss", "gradient",
                    "finite difference"), col=c("black", "red", "green"), 
                   lty=c(1,2,2))
        }
        
        
        nll <- function(y, rho=0.9){
            -lgamma((exp(y)/(exp(y) + 1)) * ((1-rho)/rho))
        }
        
        gr <- function(y, rho=0.9){
            -digamma((exp(y)/(exp(y) + 1)) * ((1-rho)/rho)) * ((1-rho)/rho) * 
                (exp(y)/(exp(y) + 1)^2)
        }
        
        
        nll <- function(y, rho=0.9){
            -lgamma((-1/(exp(y) + 1)) * ((rho-1)/rho))
        }
        
        gr <- function(y, rho=0.9){
            -digamma((-1/(exp(y) + 1)) * ((rho-1)/rho)) * ((rho-1)/rho) * 
                (exp(y)/(exp(y) + 1)^2)
        }
        
        
        x <- -100:100
        myrho <- 0.1
        par(mfrow=c(1,2))
        plot((x), nll(x, rho=myrho))
        plot((x),  gr(x, rho=myrho))
        
        x
        
        r1fun <- function(rhoi){(1 - rhoi)/rhoi}
        r2fun <- function(rhoi){(rhoi - 1)/rhoi}
        x <- seq(0,1,by=0.001)[-1][-1000]
        plot(x, r1fun(x), log="y")
        lines(x, -r2fun(x), col="red")
        
        lgamma(1e-323)
        
        #
        # Gradient and loss for E (calculated in a matrix fashion)
        #
        
        lossE <- function(par, K, N, x, wd, b, rho){
            E <- matrix(par, ncol=ncol(wd))
            y <- t(t(x %*% E %*% t(wd)) + b)
            y <- pmin(pmax(y, -700), 700)
            p <- plogis(y)
        
            rhoa  <- (1-rho)/rho
            rhob  <- (rho-1)/rho
            ey <- exp(y)
            v  <- ey/(1+ey)^2
        
        
            alpha  <- t(lgamma(t(p) * rhoa + 1e-6))
            alphaK <- t(lgamma(t(p) * rhoa + K + 0.5))
            beta   <- t(lgamma(t(p - 1) * rhob + 1e-6))
            betaNK <- t(lgamma(t(p - 1) * rhob + N - K + 0.5))
        
            mean(alpha + beta - alphaK - betaNK)
        }
        
        gradE <- function(par, x, wd, K, N, b, rho){
            E <- matrix(par, ncol=ncol(wd))
            y <- t(t(x %*% E %*% t(wd)) + b)
            y <- pmin(pmax(y, -700), 700)
            p <- plogis(y)
        
            rhoa  <- (1-rho)/rho
            rhob  <- (rho-1)/rho
            ey <- exp(y)
            v  <- ey/(1+ey)^2
        
            alpha  <- t(digamma(t(p) * rhoa + 1e-6) * rhoa)
            alphaK <- t(digamma(t(p) * rhoa + K + 0.5) * rhoa)
            beta   <- t(digamma(t(p-1) * rhob + 1e-6) * rhob)
            betaNK <- t(digamma(t(p-1) * rhob + (N - K) + 0.5) * rhob)
        
            a1 <- (v * (alpha + beta - alphaK - betaNK))
            a2 <- t(x) %*% a1
            a3 <- a2 %*% wd
            a4 <- a3 / prod(dim(y))
            grad_e <- a4
            # grad_e <- (t(x) %*%  %*% wd) / prod(dim(y))
            grad_e
        }
        
        x <- getX(K, N)
        nlgv <- numericLossGrad(lossE, 1e-8, as.vector(we), x=x, wd=wd, K=K,
                N=N, b=b, rho=rho)
        grv  <- gradE(par=as.vector(we), x=x, wd=wd, K=K, N=N, b=b, rho=rho)
        cppnlgv <- numericLossGrad(fn=truncNLL_e, epsilon=1e-8, 
                w=as.vector(we), x=x, D=wd, k=K, n=N, b=b, rho=rho[,1])
        cppgrv  <- as.vector(truncGrad_e(par=as.vector(we), x=x, D=wd, k=K, 
                n=N, b=b, rho=rho[,1]))
        
        
        plot(cppnlgv, log10(abs(cppnlgv - nlgv))); grid(); abline(0,1,col="red")
        plot(log10(sort(abs(nlgv - grv))))
        plot(log10(sort(abs(cppnlgv - cppgrv))), main=date())
        
        
    }
})
