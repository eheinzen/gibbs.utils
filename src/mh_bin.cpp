#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mh_binom(NumericVector p, NumericVector proposal, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision) {
  NumericVector out = clone(p);
  LogicalVector accept(p.size());
  for(int i=0; i < p.size(); i++) {
    double ratio = binom_LL(proposal[i], k[i], n[i], mean[i], precision[i]);
    ratio -= binom_LL(p[i], k[i], n[i], mean[i], precision[i]);
    bool a = (log(R::runif(0.0, 1.0)) <= ratio);
    accept[i] = a;
    if(accept[i]) {
      out[i] = proposal[i];
    }
  }
  out.attr("accept") = accept;
  return out;
}


// [[Rcpp::export]]
NumericVector mh_binom_mv(NumericMatrix p, NumericMatrix proposal, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q) {
  NumericMatrix out = clone(p);
  LogicalMatrix accept(p.nrow(), p.ncol());
  for(int r=0; r < p.nrow(); r++) {
    for(int i = 0; i < p.ncol(); i++) {
      double ratio = binom_LL_mv(proposal(r, _), k(r, _), n(r, _), mean(r, _), Q, i);
      ratio -= binom_LL_mv(p(r, _), k(r, _), n(r, _), mean(r, _), Q, i);
      bool a = log(R::runif(0.0, 1.0)) <= ratio;
      accept(r, i) = a;
      if(a) {
        out(r, i) = proposal(r, i);
      }
    }
  }
  out.attr("accept") = accept;
  return out;
}

