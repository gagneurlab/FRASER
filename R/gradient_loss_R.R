
#' @export
getX <- function(K, N){
    t(qlogis((K + 0.5)/(N + 1)))
}

#' @export
numericLossGrad <- function(fn, epsilon, w, ..., BPPARAM=bpparam()){
    grad <- unlist(bplapply(seq_along(w), fn=fn, w=w, epsilon=epsilon, ..., BPPARAM=BPPARAM, function(i, fn, epsilon, w, ...){
        eps <- integer(length(w))
        eps[i] <- epsilon
        (fn(w + eps, ...) - fn(w - eps, ...))/(2*epsilon)
    }))
    return(grad)
}


#
# hack to make sure its not run during package loading
#
if(FALSE){

library(BiocParallel)
library(stats)
library(VGAM)
library(devtools)

load_all()

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

N <- matrix(rnbinom(n*j, mu=cov, size=rep(disp, each=n)), ncol=n, byrow=FALSE)
K <- matrix(rbetabinom(n*j, size=N, prob=prob, rho=rho), ncol=n, byrow=FALSE)

register(SerialParam())

#
# Gradient and loss for D and b (is calculated per gene)
#
mybprange <- 1:j
res <- bplapply(mybprange, K=K, N=N, wd=wd, we=we, b=b, rho=rho, function(idx, K, N, rho, b, wd, we){
    message(idx)
    library(BiocParallel)
    # the functinos are inside to go around the windows parallelization problem (snowparam does not export globelenv :(

    maxY <- 700
    minU <- 1e-20

    ki   <- K[idx,]
    ni   <- N[idx,]
    x    <- getX(K, N)
    H    <- x %*% we
    rhoi <- rho[idx]
    bi   <- b[idx]
    wdi  <- wd[idx,]
    par <- c(bi, wdi)

    loosDb <- function(par, ki, ni, H, rhoi){
        wdi <- par[2:length(par)]
        bi <- par[1]

        yi <- H %*% wdi + bi
        yi <- pmax(pmin(yi, maxY), -maxY)
        ey <- exp(yi)

        ar <- (1 - rhoi)/rhoi
        br <- (rhoi - 1)/rhoi
        ph <- ey/(1+ey)
        u  <- -1/(ey + 1) * br
        v  <- ey/(1+ey)^2

        alpha  <- lgamma(ph * ar + minU)
        alphaK <- lgamma(ph * ar + ki + 0.5)
        beta   <- lgamma(u + minU)
        betaNK <- lgamma(u + ni - ki + 0.5)

        mean(alpha + beta - alphaK - betaNK)
    }

    gradDb <- function(par, ki, ni, H, rhoi){
        wdi <- par[2:length(par)]
        bi  <- par[1]

        yi <- H %*% wdi + bi
        yi <- pmax(pmin(yi, maxY), -maxY)
        ey <- exp(yi)

        ar <- (1 - rhoi)/rhoi
        br <- (rhoi - 1)/rhoi
        ph <- ey/(1+ey)
        v  <- ey/(1+ey)^2
        u  <- -1/(ey + 1) * br


        alpha  <- digamma(ph * ar + minU) * ar
        alphaK <- digamma(ph * ar + ki + 0.5) * ar
        beta   <- digamma(u + minU) * br
        betaNK <- digamma(u + ni - ki + 0.5) * br

        grad_b <- mean((alpha + beta - alphaK - betaNK) * v)
        grad_d <- colMeans(matrix((alpha + beta - alphaK - betaNK) * v, 100, 3) * H)

        c(grad_b, grad_d)
    }

    nlgv <- numericLossGrad(loosDb, 1e-8, par, ki=ki, ni=ni, H=H, rhoi=rhoi, BPPARAM=SerialParam())
    grv  <- gradDb(par, ki, ni, H, rhoi)
    cppnlgv <- numericLossGrad(truncNLL_db, 1e-8, par, k=ki, n=ni, H=H, rho=rhoi, BPPARAM=SerialParam())
    cppgrv  <- as.vector(truncGrad_db(par=par, H=H, k=ki, n=ni, rho=rhoi))

    list(nlgv=nlgv, grv=grv, cppnlgv=cppnlgv, cppgrv=cppgrv)
})

nlgv <- sapply(res, "[[", "nlgv")
grv  <- sapply(res, "[[", "grv")
cppnlgv <- sapply(res, "[[", "cppnlgv")
cppgrv  <- sapply(res, "[[", "cppgrv")

nlgv[,idx]
grv[,idx]

par(mfrow=c(1,2))
plot(log10(sort(abs(nlgv - grv))), ylim=c(-11,1))
plot(log10(sort(abs(cppnlgv - cppgrv))), ylim=c(-11,1))
plot(nlgv, cppnlgv); grid(); abline(0,1, col="red")

which(abs(nlgv - grv) > 0.01, arr.ind = TRUE)

nll <- function(y, rho=0.9){
    -lgamma((exp(y)/(exp(y) + 1)) * ((1-rho)/rho))
}

gr <- function(y, rho=0.9){
    -digamma((exp(y)/(exp(y) + 1)) * ((1-rho)/rho)) * ((1-rho)/rho) * (exp(y)/(exp(y) + 1)^2)
}


nll <- function(y, rho=0.9){
    -lgamma((-1/(exp(y) + 1)) * ((rho-1)/rho))
}

gr <- function(y, rho=0.9){
    -digamma((-1/(exp(y) + 1)) * ((rho-1)/rho)) * ((rho-1)/rho) * (exp(y)/(exp(y) + 1)^2)
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



library(microbenchmark)
bres <- microbenchmark(
    rl = loosDb(par, ki=ki, ni=ni, H=H, rhoi=rhoi),
    rg = gradDb(par, ki, ni, H, rhoi),
    cl = cppnll(par, k=ki, n=ni, H=H, rho=rhoi),
    cg = cppgr(par=par, H=H, k=ki, n=ni, rho=rhoi),
    times=400
)
library(ggplot2)
autoplot(bres)








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
nlgv <- numericLossGrad(lossE, 1e-8, as.vector(we), x=x, wd=wd, K=K, N=N, b=b, rho=rho, BPPARAM=MulticoreParam(10))
grv  <- gradE(par=as.vector(we), x=x, wd=wd, K=K, N=N, b=b, rho=rho)
cppnlgv <- numericLossGrad(fn=truncNLL_e, epsilon=1e-8, w=as.vector(we), x=x, D=wd, k=K, n=N, b=b, rho=rho, BPPARAM=SerialParam())
cppgrv  <- as.vector(truncGrad_e(par=as.vector(we), x=x, D=wd, k=K, n=N, b=b, rho=rho))


plot(cppnlgv, log10(abs(cppnlgv - nlgv))); grid(); abline(0,1,col="red")
plot(log10(sort(abs(nlgv - grv))))
plot(log10(sort(abs(cppnlgv - cppgrv))), main=date())








# keep the brace here to close the falls from the top
}







