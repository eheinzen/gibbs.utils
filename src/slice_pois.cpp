#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_pois_slice(double L, double k, double mean, double precision, double w, int nexpand, int ncontract) {
  double y0 = pois_LL(L, k, mean, precision) - R::rexp(1.0);
  double left = L - w;
  double right = L + w;
  int i = 0;
  while(pois_LL(left, k, mean, precision) > y0 && i++ < nexpand) {
    left -= w;
  }
  i = 0;
  while(pois_LL(right, k, mean, precision) > y0 && i++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  i = 0;
  while(pois_LL(newx, k, mean, precision) < y0 && i++ < ncontract) {
    if(newx < L) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(i == ncontract) return L;
  return newx;
}


// [[Rcpp::export]]
NumericVector slice_sample_pois(NumericVector L, NumericVector k, LogicalVector k_na, NumericVector mean, NumericVector precision,
                                double w, int nexpand, int ncontract) {
  NumericVector out(L.size());
  for(int i=0; i < L.size(); i++) {
    if(k_na[i]) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      continue;
    }
    out[i] = one_pois_slice(L[i], k[i], mean[i], precision[i], w, nexpand, ncontract);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_pois_mv(NumericMatrix L, NumericMatrix k, LogicalMatrix k_na, NumericMatrix mean, NumericMatrix Q,
                                   LogicalVector use_norm, NumericMatrix norm, double w, int nexpand, int ncontract) {
  NumericMatrix out = clone(L);
  int nrm = 0;
  for(int r=0; r < L.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      continue;
    }
    NumericVector kk = k(r, _);
    LogicalVector kk_na = k_na(r, _);
    NumericVector mm = mean(r, _);
    for(int i=0; i < L.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);
      if(kk_na[i]) {
        out(r, i) = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        continue;
      }
      out(r, i) = one_pois_slice(out(r, i), kk[i], mmm, Q(i, i), w, nexpand, ncontract);
    }
  }

  return out;
}
