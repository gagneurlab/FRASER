// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

const double epsilonU      = 1e-20;
const double MAX_EXP_VALUE = 700;

arma::mat trimmVal(arma::mat y){
    arma::uvec idx;

    idx = find(y > MAX_EXP_VALUE);
    y.elem(idx).fill(MAX_EXP_VALUE);

    idx = find(y < -MAX_EXP_VALUE);
    y.elem(idx).fill(-MAX_EXP_VALUE);

    return y;
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

    alpha  = arma::lgamma(p * rhoa + epsilonU);
    alphaK = arma::lgamma(p * rhoa + k + 0.5);
    beta   = arma::lgamma(u + epsilonU);
    betaNK = arma::lgamma(u + n - k + 0.5);

    nll = arma::accu(alpha + beta - alphaK - betaNK)/k.n_elem;

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

    alpha  = rcppdigamma(p * rhoa + epsilonU) * rhoa;
    alphaK = rcppdigamma(p * rhoa + k + 0.5) * rhoa;
    beta   = rcppdigamma(u + epsilonU) * rhob;
    betaNK = rcppdigamma(u + n - k + 0.5) * rhob;

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
    alpha  = arma::lgamma(aT + epsilonU);
    alphaK = arma::lgamma(aT + k.t() + 0.5);
    beta   = arma::lgamma(u + epsilonU);
    betaNK = arma::lgamma(u + n.t() - k.t() + 0.5);

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
    akT = rcppdigammamat(aT + k.t() + 0.5);
    aT  = rcppdigammamat(aT + epsilonU);
    bT  = rcppdigammamat(u + epsilonU);
    bnkT = rcppdigammamat(u + n.t() - k.t() + 0.5);

    alpha  = aT.each_row() % rhoa.t();
    alphaK = akT.each_row() % rhoa.t();
    beta   = bT.each_row() % rhob.t();
    betaNK = bnkT.each_row() % rhob.t();

    gr = (alpha + beta - alphaK - betaNK) % v;
    gr = (x.t() * gr * D)/y.n_elem;

    return gr;
}
