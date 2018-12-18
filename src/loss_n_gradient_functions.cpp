// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>


void printv(SEXP x){
    Rf_PrintValue(x);
}

void printmat(arma::mat x){
    printv(Rcpp::wrap(x.n_cols));
    printv(Rcpp::wrap(x.n_rows));

    x = x.head_cols(3);
    x = x.head_rows(3);
    printv(Rcpp::wrap(x));
    return;
}

void printrv(arma::rowvec v){
    printv(Rcpp::wrap(v.n_elem));
    printv(Rcpp::wrap(v.head(3)));
}

void printd(double x){
    printv(Rcpp::wrap(x));
}


const double MAX_EXP_VALUE = 700;

arma::mat minValForExp(arma::mat y){
    arma::uvec idx;

    idx = find(y > MAX_EXP_VALUE);
    y.elem(idx).fill(MAX_EXP_VALUE);

    idx = find(y < -MAX_EXP_VALUE);
    y.elem(idx).fill(-MAX_EXP_VALUE);

    return y;
}

// [[Rcpp::export()]]
arma::mat predictMatY(arma::mat x, arma::mat E, arma::mat D, arma::vec b){
    arma::mat y = x * E * D.t();
    y.each_row() += b.t();
    y = minValForExp(y);

    return y;
}

// [[Rcpp::export()]]
arma::mat predictMatC(arma::mat x, arma::mat E, arma::mat D, arma::vec b,
                    arma::vec sf){
    arma::mat y = predictMatY(x, E, D, b) ;
    arma::mat c = arma::exp(y);
    c.each_col() %= sf;

    return c;
}

arma::vec colMeans(arma::mat X){
    return arma::vectorise(arma::sum(X,0))/X.n_rows;
}

// [[Rcpp::export()]]
double truncNLL_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho){
    double b, rhoa;
    arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, ll;

    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);

    y = H * d + b;
    y = minValForExp(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    p    = ey/(1 + ey);
    u    = -1/(1 + ey);
    v    = ey/arma::sqrt(1 + ey);

    alpha  <- arma::lgamma(p * rhoa);
    alphaK <- arma::lgamma(p * rhoa + k + 0.5);
    beta   <- arma::lgamma(u);
    betaNK <- arma::lgamma(u + n - k + 0.5);

    ll = arma::accu(alpha + beta - alphaK - betaNK)/k.n_elem;

    printmat(y);
    printd(b);

    return arma::as_scalar(-ll);
}

// [[Rcpp::export()]]
arma::vec truncGrad_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho){
    double b, rhoa, rhob;
    arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, grb, grd;

    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);

    y = H * d + b;
    y = minValForExp(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(1 + ey);
    v    = ey/arma::sqrt(1 + ey);

    alpha  <- arma::lgamma(p * rhoa) * rhoa;
    alphaK <- arma::lgamma(p * rhoa + k + 0.5) * rhoa;
    beta   <- arma::lgamma(u) * rhob;
    betaNK <- arma::lgamma(u + n - k + 0.5) * rhob;

    grb = arma::accu((alpha + beta - alphaK - betaNK) % v)/k.n_elem;
    grd = (alpha + beta - alphaK - betaNK) % v;
    grd = colMeans(grd % H.each_col());

    arma::mat ans = arma::join_cols(grb, grd);
    return ans.col(0);
}
