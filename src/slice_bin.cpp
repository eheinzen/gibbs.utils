#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_binom_slice(double p, double k, double n, double mean, double precision, double w, int nexpand, int ncontract) {
  double y0 = binom_LL(p, k, n, mean, precision) - R::rexp(1.0);
  double left = p - w;
  double right = p + w;
  int i = 0;
  while(binom_LL(left, k, n, mean, precision) > y0 && i++ < ncontract) {
    left -= w;
  }
  i = 0;
  while(binom_LL(right, k, n, mean, precision) > y0 && i++ < ncontract) {
    right += w;
  }
  double newx = R::runif(left, right);
  i = 0;
  while(binom_LL(newx, k, n, mean, precision) < y0 && i++ < ncontract) {
    if(newx < p) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(i == ncontract) return p;
  return newx;
}




// [[Rcpp::export]]
NumericVector slice_sample_binom(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision,
                                 double w, int nexpand, int ncontract) {
  NumericVector out(p.size());
  for(int i=0; i < p.size(); i++) {
    if(n[i] == 0.0) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      continue;
    }
    out[i] = one_binom_slice(p[i], k[i], n[i], mean[i], precision[i], w, nexpand, ncontract);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix slice_sample_binom_mv(NumericMatrix p, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q,
                                    LogicalVector use_norm, NumericMatrix norm, double w, int nexpand, int ncontract) {
  NumericMatrix out = clone(p);
  int nrm = 0;
  for(int r=0; r < p.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      continue;
    }
    NumericVector kk = k(r, _);
    NumericVector nn = n(r, _);
    NumericVector mm = mean(r, _);
    for(int i=0; i < p.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);
      if(nn[i] == 0.0) {
        out(r, i) = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        continue;
      }
      out(r, i) = one_binom_slice(out(r, i), kk[i], nn[i], mmm, Q(i, i), w, nexpand, ncontract);
    }
  }
  return out;
}
