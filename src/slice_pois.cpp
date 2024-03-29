#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_pois_slice(const double L, const double k,
                      const double mean, const double precision,
                      const double trunc_at, const bool lower,
                      const double w, const int nexpand, const int ncontract) {
  double y0 = pois_LL(L, k, mean, precision, trunc_at, lower) - R::rexp(1.0);
  double left = L - w;
  double right = L + w;
  int i = 0;
  while(pois_LL(left, k, mean, precision, trunc_at, lower) > y0 && i++ < nexpand) {
    left -= w;
  }
  i = 0;
  while(pois_LL(right, k, mean, precision, trunc_at, lower) > y0 && i++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  i = 0;
  while(pois_LL(newx, k, mean, precision, trunc_at, lower) < y0 && i++ < ncontract) {
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
NumericVector slice_sample_pois(const NumericVector L,
                                const NumericVector k, const LogicalVector k_na,
                                const NumericVector mean, const LogicalVector mean_inf, const NumericVector precision,
                                const NumericVector trunc_at, const LogicalVector lower,
                                const double w, const int nexpand, const int ncontract) {
  NumericVector out(L.size());
  for(int i=0; i < L.size(); i++) {
    if(k_na[i]) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      continue;
    } else if(mean_inf[i]) {
      out[i] = R_NegInf;
      continue;
    }
    out[i] = one_pois_slice(L[i], k[i], mean[i], precision[i], trunc_at[i], lower[i], w, nexpand, ncontract);
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_pois_mv(const NumericMatrix L,
                                   const NumericMatrix k, const LogicalMatrix k_na,
                                   const NumericMatrix mean, const NumericMatrix Q,
                                   const NumericMatrix trunc_at, const LogicalMatrix lower,
                                   const LogicalVector use_norm, const NumericMatrix norm,
                                   const double w, const int nexpand, const int ncontract) {
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
    NumericVector aa = trunc_at(r, _);
    LogicalVector bb = lower(r, _);
    NumericVector tmp = out(r, _);
    for(int i=0; i < L.ncol(); i++) {
      // safe not to replace tmp[i] because it's not used in this calculation
      double mmm = cond_mv_mean(tmp, mm, Q, i);
      if(kk_na[i]) {
        tmp[i] = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        continue;
      }
      tmp[i] = one_pois_slice(tmp[i], kk[i], mmm, Q(i, i), aa[i], bb[i], w, nexpand, ncontract);
    }
    out(r, _) = tmp;
  }

  return out;
}
