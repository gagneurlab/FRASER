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
    if(x.n_rows >= 1){
        x = x.head_rows(5);
        x = x.head_cols(5);
        printv(Rcpp::wrap(x));
    } else {
        prints("No rows exists here!");
    }
    return;
}

void printrv(arma::rowvec v){
    printv(Rcpp::wrap(v));
}

void printvec(arma::vec v){
    printv(Rcpp::wrap(v));
    //.head(std::min(v.n_elem-0.1, 5 - 0.1) + 0.1))
}

void printd(double x){
    printv(Rcpp::wrap(x));
}
