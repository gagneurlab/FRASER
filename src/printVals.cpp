#include <RcppArmadillo.h>

using namespace Rcpp;

const bool printVals = true;

void printv(SEXP x){
    if(printVals == true){
        print(x);
    }
}

void prints(std::string s){
    printv(Rcpp::wrap(s));
}

void printmat(arma::mat x){
    printv(Rcpp::wrap(x.n_cols));
    printv(Rcpp::wrap(x.n_rows));
    if((x.n_rows >= 1) & (x.n_cols >= 1)){
        int nCols = as<int>(wrap(x.n_cols));
        int nRows = as<int>(wrap(x.n_rows));
        nCols = std::min<int>(nCols, 5);
        nRows = std::min<int>(nRows, 5);

        x = x.head_rows(nRows);
        x = x.head_cols(nCols);
        printv(Rcpp::wrap(x));
    } else {
        prints("No rows exists here!");
    }
    return;
}

void printvec(arma::vec v){
    int nElem = as<int>(wrap(v.n_elem));
    nElem = std::min<int>(nElem, 20);
    printv(Rcpp::wrap(v.head(nElem)));
}

void printrv(arma::rowvec v){
    printvec(arma::vectorise(v));
}

void printd(double x){
    printv(Rcpp::wrap(x));
}
