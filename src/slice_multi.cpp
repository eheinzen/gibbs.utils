#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_multinom_slice(NumericVector p_j, LogicalVector z_j, double k, double n, double mean, double precision, int j, double w, int nexpand, int ncontract) {
  if(n == 0.0 || !z_j[j]) {
    return R::rnorm(mean, 1/sqrt(precision));
  }
  double y0 = multinom_LL(p_j, z_j, k, n, mean, precision, j) - R::rexp(1.0);
  double pj = p_j[j];
  double left = pj - w;
  double right = pj + w;
  int jj = 0;
  while(multinom_LL(replace_it(p_j, j, left), z_j, k, n, mean, precision, j) > y0 && jj++ < nexpand) {
    left -= w;
  }
  jj = 0;
  while(multinom_LL(replace_it(p_j, j, right), z_j, k, n, mean, precision, j) > y0 && jj++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  jj = 0;
  while(multinom_LL(replace_it(p_j, j, newx), z_j, k, n, mean, precision, j) < y0 && jj++ < ncontract) {
    if(newx < pj) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(jj == ncontract) return pj;
  return newx;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_multinom_mv(List p_j, LogicalMatrix z, NumericMatrix k, NumericMatrix n, NumericMatrix p_i, NumericMatrix mean, NumericMatrix Q, int j, double w, int nexpand, int ncontract) {
  // p_j is a list of length r with matrices of dim i x j
  // p_i is a matrix of dim r x i
  // because we don't loop through the j-dim here, and because each r is assumed independent,
  //   we don't reuse and elements of p_j in the log-likelihood; therefore, no need to update in this step
  NumericMatrix out = clone(p_i);
  for(int r=0; r < p_i.nrow(); r++) {
    NumericMatrix pp_jj = p_j[r];
    NumericVector kk = k(r, _);
    NumericVector nn = n(r, _);
    NumericVector mm = mean(r, _);
    for(int i=0; i < p_i.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);
      out(r, i) = one_multinom_slice(pp_jj(i, _), z(i, _), kk[i], nn[i], mmm, Q(i, i), j, w, nexpand, ncontract);
    }
  }
  return out;
}
