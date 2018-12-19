// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "printVals.h"

using namespace Rcpp;


const double epsilonU      = 0.000001;
const double MAX_EXP_VALUE = 700;

const double MAX_DIG_VALUE = 36.7369;

arma::mat trimmVal(arma::mat y, double val=MAX_EXP_VALUE,
                    bool neg=true, bool pos=true){
    arma::uvec idx;

    if(pos == true){
        idx = find(y > val);
        y.elem(idx).fill(val);
    }

    if(neg == true){
        idx = find(y < -val);
        y.elem(idx).fill(-val);
    }

    return y;
}


arma::mat replaceInf(arma::mat x, double val=MAX_DIG_VALUE){
    arma::uvec idx;


    idx = arma::find_nonfinite(x);
    x.elem(idx).fill(val);

    return x;
}


arma::vec colMeans(arma::mat X){
    return arma::vectorise(arma::sum(X,0))/X.n_rows;
}


arma::vec rcppdigamma(arma::vec x){
    NumericVector out = as<NumericVector>(wrap(x));
    out = digamma(out);
    return as<arma::vec>(wrap(out));
}

arma::mat rcppdigammamat(arma::mat x){
    NumericVector out = as<NumericVector>(wrap(arma::vectorise(x)));
    out = digamma(out);

    arma::vec myvec = as<arma::vec>(wrap(out));
    arma::mat outmat = arma::reshape(myvec, x.n_rows, x.n_cols);
    return outmat;
}


// [[Rcpp::export()]]
double truncNLL_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho){
    double b, rhoa, rhob;
    arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll;

    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);

    y = H * d + b;
    y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(1 + ey) * rhob;
    v    = ey/arma::square(1 + ey);

    alpha  = arma::lgamma(p * rhoa);
    alphaK = arma::lgamma(p * rhoa + k);
    beta   = arma::lgamma(u + epsilonU);
    betaNK = arma::lgamma(u + n - k + 0.5);

    beta   <- replaceInf(beta);
    betaNK <- replaceInf(betaNK);

    nll = arma::accu(alpha + beta - alphaK - betaNK)/k.n_elem;

    //printd(b);
    //printvec(d);
    //printmat(ey);
    //printmat(v);
    //printmat(p);
    //prints("alpha");
    //printvec(u);
    //printvec(beta);
    //printvec(betaNK);

    return arma::as_scalar(nll);
}


// [[Rcpp::export()]]
arma::vec truncGrad_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho){
    double b, rhoa, rhob;
    arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, grb, grd;

    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);

    y = H * d + b;
    y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(1 + ey) * rhob;
    v    = ey/arma::square(1 + ey);

    alpha  = rcppdigamma(p * rhoa) * rhoa;
    alphaK = rcppdigamma(p * rhoa + k) * rhoa;
    beta   = rcppdigamma(u + epsilonU) * rhob;
    betaNK = rcppdigamma(u + n - k + 0.5) * rhob;

    beta <- replaceInf(beta);
    betaNK <- replaceInf(betaNK);

    grb = arma::accu((alpha + beta - alphaK - betaNK) % v)/k.n_elem;
    grd = (alpha + beta - alphaK - betaNK) % v;
    grd = colMeans(grd % H.each_col());

    arma::mat ans = arma::join_cols(grb, grd);
    return ans.col(0);
}

// [[Rcpp::export()]]
double truncNLL_e(arma::vec par, arma::mat x, arma::mat D, arma::vec b,
                    arma::mat k, arma::mat n, arma::vec rho){
    arma::vec rhoa, rhob;
    arma::mat E, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll, aT;

    E = arma::reshape(par, D.n_rows, D.n_cols);

    y = x * E * D.t();
    y = y.each_row() + b.t();
    y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(ey + 1);
    u    = u.each_row() % rhob.t();
    v    = ey/arma::square(1 + ey);

    aT = p.each_row() % rhoa.t();
    alpha  = arma::lgamma(aT);
    alphaK = arma::lgamma(aT + k.t());
    beta   = arma::lgamma(u + epsilonU);
    betaNK = arma::lgamma(u + n.t() - k.t() + 0.5);

    beta   <- replaceInf(beta);
    betaNK <- replaceInf(betaNK);

    nll = arma::accu(alpha + beta - alphaK - betaNK)/k.n_elem;

    return arma::as_scalar(nll);
}


// [[Rcpp::export()]]
arma::mat truncGrad_e(arma::vec par, arma::mat x, arma::mat D, arma::vec b,
                    arma::mat k, arma::mat n, arma::vec rho){
    arma::vec rhoa, rhob;
    arma::mat E, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll, aT, akT, bT, bnkT, gr;

    E = arma::reshape(par, D.n_rows, D.n_cols);

    y = x * E * D.t();
    y = y.each_row() + b.t();
    y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(ey + 1);
    u    = u.each_row() % rhob.t();
    v    = ey/arma::square(1 + ey);

    aT  = p.each_row() % rhoa.t();
    akT = rcppdigammamat(aT + k.t());
    aT  = rcppdigammamat(aT);
    bT  = rcppdigammamat(u + epsilonU);
    bnkT = rcppdigammamat(u + n.t() - k.t() + 0.5);

    alpha  = aT.each_row() % rhoa.t();
    alphaK = akT.each_row() % rhoa.t();
    beta   = bT.each_row() % rhob.t();
    betaNK = bnkT.each_row() % rhob.t();

    beta <- replaceInf(beta);
    betaNK <- replaceInf(betaNK);

    gr = (alpha + beta - alphaK - betaNK) % v;
    gr = (x.t() * gr * D);
    gr = gr / (D.n_cols * D.n_rows);

    return gr;
}
