// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

const double MAX_EXP_VALUE = 700;
double PSEUDO_COUNT = 1;

void setPseudoCount(double pseudoCount){
    if(pseudoCount >= 0){
        PSEUDO_COUNT = pseudoCount;
    }
}

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

arma::vec rowMeans(arma::mat X){
    return arma::vectorise(arma::sum(X,1))/X.n_cols;
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

arma::vec estMu(arma::vec y, arma::uvec pos){
  arma::vec ypos, est;
  ypos = y.elem( pos);
  est.ones(pos.n_elem);
  est.elem( arma::find(ypos < 0) ).zeros();
  return est;
}

// [[Rcpp::export()]]
arma::mat predictYCpp(arma::mat H, arma::mat D, arma::vec b){

  arma::mat y, ey, mu;

  y = H * D.t();
  y = y.each_row() + b.t();

  return y;
}

// [[Rcpp::export()]]
arma::mat predictMuCpp(arma::mat y){

  arma::mat ey, mu;

  ey = arma::exp(y);
  mu = ey/(1 + ey);

  mu.elem( arma::find_nonfinite(mu) ) = estMu(vectorise(y), arma::find_nonfinite(mu));

  return mu;
}

// [[Rcpp::export()]]
arma::vec estLgammaAlpha(arma::vec y, arma::uvec pos, double ar){
  arma::vec ypos, est;
  ypos = y.elem( pos);
  arma::mat mu(1,1);
  mu = -35.0;
  double fixed = arma::as_scalar(lgamma(predictMuCpp(mu) * ar));
  est = fixed - 1*(ypos+35);
  return est;
}

// [[Rcpp::export()]]
arma::vec estLgammaBeta(arma::vec y, arma::uvec pos, double br){
  arma::vec ypos, est;
  ypos = y.elem( pos);
  arma::mat mu(1,1);
  mu = 30.0;
  double fixed = arma::as_scalar(lgamma((predictMuCpp(mu)-1) * br));
  est = fixed + 1*(ypos-30);
  return est;
}

// [[Rcpp::export()]]
double truncNLL_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho, double lambda){
  double b, rhoa, rhob;
  arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll;

  b = par.at(0);
  d = par.subvec(1, par.n_elem-1);

  y = H * d + b;
  //y = trimmVal(y);
  ey = arma::exp(y);

  rhoa = (1 - rho)/rho;
  rhob = (rho - 1)/rho;
  p    = ey/(1 + ey);
  u    = -1/(1 + ey);

  p.elem( arma::find_nonfinite(p) ) = estMu(y, arma::find_nonfinite(p));
  arma::uvec uPos = arma::find_nonfinite(u);
  u.elem( uPos ) = arma::repelem(p-1, 1, uPos.n_elem );

  u    = u * rhob;

  alpha  = arma::lgamma(p * rhoa );
  alphaK = arma::lgamma(p * rhoa + k + PSEUDO_COUNT);
  beta   = arma::lgamma(u);
  betaNK = arma::lgamma(u + n - k + PSEUDO_COUNT);

  // arma::vec abs;
  arma::uvec infPosA, infPosB;
  // abs = arma::abs(y);
  infPosA = arma::find_nonfinite(alpha);
  // alpha.elem( infPosA ) = abs.elem( infPosA );
  alpha.elem( infPosA ) = estLgammaAlpha(y, infPosA, rhoa);
  infPosB = arma::find_nonfinite(beta);
  // beta.elem( infPosB ) = abs.elem( infPosB );
  beta.elem( infPosB ) = estLgammaBeta(y, infPosB, rhob);

  nll = arma::accu(alpha + beta - alphaK - betaNK)/k.n_elem;

  nll = nll + (lambda/k.n_elem) * arma::accu(d % d);

  return arma::as_scalar(nll);
}

arma::vec estDigammaAlpha(arma::vec y, arma::uvec pos){
  arma::vec ypos, est, ones;
  ypos = y.elem( pos);
  est.zeros(pos.n_elem);
  arma::uvec yNegPos;
  yNegPos = arma::find(ypos < 0);
  est.elem( yNegPos ) = -ones.ones(yNegPos.n_elem); // = -1
  return est;
}

arma::vec estDigammaBeta(arma::vec y, arma::uvec pos){
  arma::vec ypos, est;
  ypos = y.elem( pos);
  est.ones(pos.n_elem);
  est.elem( arma::find(ypos < 0) ).zeros();
  return est;
}

// [[Rcpp::export()]]
arma::vec truncGrad_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho, double lambda){
    double b, rhoa, rhob;
    arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, grb, grd;

    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);

    y = H * d + b;
    //y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(1 + ey);

    p.elem( arma::find_nonfinite(p) ) = estMu(y, arma::find_nonfinite(p));
    arma::uvec uPos = arma::find_nonfinite(u);
    u.elem( uPos ) = arma::repelem(p-1, 1, uPos.n_elem );

    u    = u * rhob;
    v    = ey/arma::square(1 + ey);

    alpha  = (rcppdigamma(p * rhoa) * rhoa) % v;
    alphaK = (rcppdigamma(p * rhoa + k + PSEUDO_COUNT) * rhoa) % v;
    beta   = (rcppdigamma(u) * rhob) % v;
    betaNK = (rcppdigamma(u + n - k + PSEUDO_COUNT) * rhob) % v;

    v.elem( arma::find_nonfinite(v) ).zeros();
    alpha.elem( arma::find(v == 0) ) = estDigammaAlpha(y, arma::find(v == 0));
    beta.elem( arma::find(v == 0) ) = estDigammaBeta(y, arma::find(v == 0));

    arma::uvec infPosA, infPosB, infPosAk, infPosBnk, ypos;
    infPosA = arma::find_nonfinite(alpha);
    alpha.elem( infPosA ) = estDigammaAlpha(y, infPosA);
    infPosB = arma::find_nonfinite(beta);
    beta.elem( infPosB ) = estDigammaBeta(y, infPosB);
    infPosAk = arma::find_nonfinite(alphaK);
    alphaK.elem( infPosAk ).zeros();
    infPosBnk = arma::find_nonfinite(betaNK);
    betaNK.elem( infPosBnk ).zeros();

    // grb = arma::accu((alpha + beta - alphaK - betaNK) % v)/k.n_elem;
    // grd = (alpha + beta - alphaK - betaNK) % v;
    // grd = colMeans(grd % H.each_col());

    grb = arma::accu((alpha + beta - alphaK - betaNK))/k.n_elem;
    grd = (alpha + beta - alphaK - betaNK) ;
    grd = colMeans(grd % H.each_col());

    grd = grd + (2*lambda/k.n_elem) * d;

    arma::mat ans = arma::join_cols(grb, grd);
    return ans.col(0);
}

// [[Rcpp::export()]]
double truncNLL_e(arma::vec par, arma::mat x, arma::mat D, arma::vec b,
                    arma::mat k, arma::mat n, arma::vec rho){
    arma::vec rhoa, rhob;
    arma::mat E, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll, aT, bT;

    E = arma::reshape(par, x.n_cols, D.n_cols);

    y = x * E * D.t();
    y = y.each_row() + b.t();
    //y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(ey + 1);

    p.elem( arma::find_nonfinite(p) ) = estMu(vectorise(y), arma::find_nonfinite(p));
    arma::uvec uPos = arma::find_nonfinite(u);
    u.elem( uPos ) = arma::repelem(p-1, 1, uPos.n_elem );

    v    = ey/arma::square(1 + ey);
    aT = p.each_row() % rhoa.t();
    bT = u.each_row() % rhob.t();

    alpha  = arma::lgamma(aT);
    alphaK = arma::lgamma(aT + k.t() + PSEUDO_COUNT);
    beta   = arma::lgamma(bT);
    betaNK = arma::lgamma(bT + n.t() - k.t() + PSEUDO_COUNT);

    arma::mat abs;
    arma::uvec infPosA, infPosB;
    abs = arma::abs(y);
    infPosA = arma::find_nonfinite(alpha);
    alpha.elem( infPosA ) = abs.elem( infPosA );
    // alpha.elem( infPosA ) = estLgammaAlpha(y, infPosA, rhoa);
    infPosB = arma::find_nonfinite(beta);
    beta.elem( infPosB ) = abs.elem( infPosB );
    // beta.elem( infPosB ) = estLgammaBeta(y, infPosB, rhob);

    nll = arma::accu(alpha + beta - alphaK - betaNK)/k.n_elem;

    return arma::as_scalar(nll);
}

// [[Rcpp::export()]]
arma::mat truncGrad_e(arma::vec par, arma::mat x, arma::mat D, arma::vec b,
                    arma::mat k, arma::mat n, arma::vec rho){
    arma::vec rhoa, rhob;
    arma::mat E, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll, aT, akT, bT, bnkT, gr;

    E = arma::reshape(par, x.n_cols, D.n_cols);

    y = x * E * D.t();
    y = y.each_row() + b.t();
    //y = trimmVal(y);
    ey = arma::exp(y);

    rhoa = (1 - rho)/rho;
    rhob = (rho - 1)/rho;
    p    = ey/(1 + ey);
    u    = -1/(ey + 1);

    p.elem( arma::find_nonfinite(p) ) = estMu(vectorise(y), arma::find_nonfinite(p));
    arma::uvec uPos = arma::find_nonfinite(u);
    u.elem( uPos ) = arma::repelem(p-1, 1, uPos.n_elem );

    u   = u.each_row() % rhob.t();
    v    = ey/arma::square(1 + ey);
    aT  = p.each_row() % rhoa.t();

    akT = rcppdigammamat(aT + k.t() + PSEUDO_COUNT);
    aT  = rcppdigammamat(aT);
    bT  = rcppdigammamat(u);
    bnkT = rcppdigammamat(u + n.t() - k.t() + PSEUDO_COUNT);

    // alpha  = aT.each_row() % rhoa.t();
    // alphaK = akT.each_row() % rhoa.t();
    // beta   = bT.each_row() % rhob.t();
    // betaNK = bnkT.each_row() % rhob.t();
    //
    // gr = (alpha + beta - alphaK - betaNK) % v;
    // gr = (x.t() * gr * D)/y.n_elem;

    alpha  = (aT.each_row() % rhoa.t() ) % v;
    alphaK = (akT.each_row() % rhoa.t() ) % v;
    beta   = (bT.each_row() % rhob.t() ) % v;
    betaNK = (bnkT.each_row() % rhob.t() ) % v;

    v.elem( arma::find_nonfinite(v) ).zeros();
    alpha.elem( arma::find(v == 0) ) = estDigammaAlpha(vectorise(y), arma::find(v == 0));
    beta.elem( arma::find(v == 0) ) = estDigammaBeta(vectorise(y), arma::find(v == 0));

    arma::uvec infPosA, infPosB, infPosAk, infPosBnk, ypos;
    infPosA = arma::find_nonfinite(alpha);
    alpha.elem( infPosA ) = estDigammaAlpha(vectorise(y), infPosA);
    infPosB = arma::find_nonfinite(beta);
    beta.elem( infPosB ) = estDigammaBeta(vectorise(y), infPosB);
    infPosAk = arma::find_nonfinite(alphaK);
    alphaK.elem( infPosAk ).zeros();
    infPosBnk = arma::find_nonfinite(betaNK);
    betaNK.elem( infPosBnk ).zeros();

    gr = (alpha + beta - alphaK - betaNK);
    gr = (x.t() * gr * D)/y.n_elem;

    return gr;
}

// [[Rcpp::export()]]
double truncNLL_rho(double rho, arma::vec yi, arma::vec ki, arma::vec ni){
  arma::vec mui, u, alpha, alphaK, beta, betaNK, alphaBeta, nll;
  double rhoa, rhob;

  rhoa = (1 - rho)/rho;
  rhob = (rho - 1)/rho;
  mui = predictMuCpp(yi);
  u    = (mui-1) * rhob;

  alpha  = arma::lgamma(mui * rhoa);
  alphaK = arma::lgamma(mui * rhoa + ki + PSEUDO_COUNT);
  beta   = arma::lgamma(u);
  betaNK = arma::lgamma(u + ni - ki + PSEUDO_COUNT);
  alphaBeta = arma::lgamma(rhoa + ni + (2*PSEUDO_COUNT)) - lgamma(rhoa);

  // arma::vec abs;
  arma::uvec infPosA, infPosB;
  // abs = arma::abs(yi);
  infPosA = arma::find_nonfinite(alpha);
  // alpha.elem( infPosA ) = abs.elem( infPosA );
  alpha.elem( infPosA ) = estLgammaAlpha(yi, infPosA, rhoa);
  infPosB = arma::find_nonfinite(beta);
  // beta.elem( infPosB ) = abs.elem( infPosB );
  beta.elem( infPosB ) = estLgammaBeta(yi, infPosB, rhob);

  nll = arma::accu(alpha + beta - alphaK - betaNK + alphaBeta)/ki.n_elem;

  return arma::as_scalar(nll);
}

// [[Rcpp::export()]]
arma::vec fullNLL(arma::mat y, arma::mat rho, arma::mat k, arma::mat n, arma::mat D, double lambda, bool byRows=false){
  arma::mat rhoa, rhob;
  arma::mat mu, u, alpha, alphaK, beta, betaNK, nonTruncTerms, nll, aT;


  rhoa = (1 - rho)/rho;
  rhob = -rhoa;
  mu = predictMuCpp(y);
  aT = mu % rhoa;
  u    = mu - 1;
  u    = u % rhob;

  alpha  = arma::lgamma(aT);
  alphaK = arma::lgamma(aT + k + PSEUDO_COUNT);
  beta   = arma::lgamma(u);
  betaNK = arma::lgamma(u + n - k + PSEUDO_COUNT);
  nonTruncTerms = - arma::lgamma(n+1+(2*PSEUDO_COUNT)) + arma::lgamma(k+1+PSEUDO_COUNT) + arma::lgamma(n-k+1+PSEUDO_COUNT)
                  + arma::lgamma(rhoa + n + (2*PSEUDO_COUNT)) - arma::lgamma(rhoa);

  arma::mat abs;
  arma::uvec infPosA, infPosB;
  abs = arma::abs(y);
  infPosA = arma::find_nonfinite(alpha);
  alpha.elem( infPosA ) = abs.elem( infPosA );
  // alpha.elem( infPosA ) = estLgammaAlpha(y, infPosA, rhoa);
  infPosB = arma::find_nonfinite(beta);
  beta.elem( infPosB ) = abs.elem( infPosB );
  // beta.elem( infPosB ) = estLgammaBeta(y, infPosB, rhob);

  if(byRows){
    nll = rowMeans(alpha + beta - alphaK - betaNK + nonTruncTerms);
    nll = nll + (lambda/k.n_elem) * rowMeans(D % D);
    return arma::vectorise(nll);
  } else {
    nll = arma::accu(alpha + beta - alphaK - betaNK + nonTruncTerms)/k.n_elem;
    nll = nll + (lambda/k.n_elem) * arma::accu(D % D);
    return arma::vectorise(nll);
  }
}

// get weights
arma::vec getWeights(arma::vec k, arma::vec n, arma::vec mu, double rho){
  double c;
  arma::vec r, w;
  
  // pearson residuals for BB
  r = ((k+PSEUDO_COUNT) - (n+2*PSEUDO_COUNT) % mu) / sqrt((n+2*PSEUDO_COUNT) % mu % (1-mu) % (1+((n+2*PSEUDO_COUNT)-1)*rho));
  
  // weights according to Huber function
  c = 1.345; // constant, as in edgeR
  w.ones(r.n_elem);
  arma::uvec pos = arma::find(abs(r) > c);
  w.elem(pos) = c/abs(r.elem(pos));

  return w;
  
}


// weighted NLL
// [[Rcpp::export()]]
double truncWeightedNLL_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho, double lambda, arma::vec w){
  double b, rhoa, rhob;
  arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, nll;
  
  b = par.at(0);
  d = par.subvec(1, par.n_elem-1);
  
  y = H * d + b;
  //y = trimmVal(y);
  ey = arma::exp(y);
  
  // arma::vec w;
  // w.ones(ey.n_elem);
  // if(weighted){
  //   w = getWeights(k, n, ey, rho);
  // }
  
  rhoa = (1 - rho)/rho;
  rhob = (rho - 1)/rho;
  p    = ey/(1 + ey);
  u    = -1/(1 + ey);
  
  p.elem( arma::find_nonfinite(p) ) = estMu(y, arma::find_nonfinite(p));
  arma::uvec uPos = arma::find_nonfinite(u);
  u.elem( uPos ) = arma::repelem(p-1, 1, uPos.n_elem );
  
  u    = u * rhob;
  
  alpha  = arma::lgamma(p * rhoa );
  alphaK = arma::lgamma(p * rhoa + k + PSEUDO_COUNT);
  beta   = arma::lgamma(u);
  betaNK = arma::lgamma(u + n - k + PSEUDO_COUNT);
  
  // arma::vec abs;
  arma::uvec infPosA, infPosB;
  // abs = arma::abs(y);
  infPosA = arma::find_nonfinite(alpha);
  // alpha.elem( infPosA ) = abs.elem( infPosA );
  alpha.elem( infPosA ) = estLgammaAlpha(y, infPosA, rhoa);
  infPosB = arma::find_nonfinite(beta);
  // beta.elem( infPosB ) = abs.elem( infPosB );
  beta.elem( infPosB ) = estLgammaBeta(y, infPosB, rhob);
  
  nll = arma::accu((alpha + beta - alphaK - betaNK)%w)/k.n_elem;
  
  nll = nll + (lambda/k.n_elem) * arma::accu(d % d);
  
  return arma::as_scalar(nll);
}

// weighted gradient of NLL
// [[Rcpp::export()]]
arma::vec truncWeightedGrad_db(arma::vec par, arma::mat H, arma::vec k, arma::vec n, double rho, double lambda, arma::vec w){
  double b, rhoa, rhob;
  arma::vec d, y, ey, p, u, v, alpha, alphaK, beta, betaNK, grb, grd;
  
  b = par.at(0);
  d = par.subvec(1, par.n_elem-1);
  
  y = H * d + b;
  //y = trimmVal(y);
  ey = arma::exp(y);
  
  // arma::vec w;
  // w.ones(ey.n_elem);
  // if(weighted){
  //   w = getWeights(k, n, ey, rho);
  // }
  
  rhoa = (1 - rho)/rho;
  rhob = (rho - 1)/rho;
  p    = ey/(1 + ey);
  u    = -1/(1 + ey);
  
  p.elem( arma::find_nonfinite(p) ) = estMu(y, arma::find_nonfinite(p));
  arma::uvec uPos = arma::find_nonfinite(u);
  u.elem( uPos ) = arma::repelem(p-1, 1, uPos.n_elem );
  
  u    = u * rhob;
  v    = ey/arma::square(1 + ey);
  
  alpha  = (rcppdigamma(p * rhoa) * rhoa) % v;
  alphaK = (rcppdigamma(p * rhoa + k + PSEUDO_COUNT) * rhoa) % v;
  beta   = (rcppdigamma(u) * rhob) % v;
  betaNK = (rcppdigamma(u + n - k + PSEUDO_COUNT) * rhob) % v;
  
  v.elem( arma::find_nonfinite(v) ).zeros();
  alpha.elem( arma::find(v == 0) ) = estDigammaAlpha(y, arma::find(v == 0));
  beta.elem( arma::find(v == 0) ) = estDigammaBeta(y, arma::find(v == 0));
  
  arma::uvec infPosA, infPosB, infPosAk, infPosBnk, ypos;
  infPosA = arma::find_nonfinite(alpha);
  alpha.elem( infPosA ) = estDigammaAlpha(y, infPosA);
  infPosB = arma::find_nonfinite(beta);
  beta.elem( infPosB ) = estDigammaBeta(y, infPosB);
  infPosAk = arma::find_nonfinite(alphaK);
  alphaK.elem( infPosAk ).zeros();
  infPosBnk = arma::find_nonfinite(betaNK);
  betaNK.elem( infPosBnk ).zeros();
  
  grb = arma::accu((alpha + beta - alphaK - betaNK) % w)/k.n_elem;
  grd = (alpha + beta - alphaK - betaNK) % w;
  grd = colMeans(grd % H.each_col());
  
  grd = grd + (2*lambda/k.n_elem) * d;
  
  arma::mat ans = arma::join_cols(grb, grd);
  return ans.col(0);
}