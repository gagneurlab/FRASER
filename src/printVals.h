#ifndef PrintVals_H
#define PrintVals_H

#include <RcppArmadillo.h>

using namespace Rcpp;


const bool printVals = false;

void printv(SEXP x);
void prints(std::string s);
void printmat(arma::mat x);
void printrv(arma::rowvec v);
void printvec(arma::vec v);
void printd(double x);

#endif
