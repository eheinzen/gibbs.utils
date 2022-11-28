#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_multinom_slice(NumericVector p_ij, LogicalVector z_ij, IntegerVector which_i, double k, double n, double mean, double precision, int ij, double w, int nexpand, int ncontract) {
  if(n == 0.0 || !z_ij[ij]) {
    return R::rnorm(mean, 1/sqrt(precision));
  }
  double y0 = multinom_LL(p_ij, z_ij, which_i, k, n, mean, precision, ij) - R::rexp(1.0);
  double pij = p_ij[ij];
  double left = pij - w;
  double right = pij + w;
  int jj = 0;
  while(multinom_LL(replace_it(p_ij, ij, left), z_ij, which_i, k, n, mean, precision, ij) > y0 && jj++ < nexpand) {
    left -= w;
  }
  jj = 0;
  while(multinom_LL(replace_it(p_ij, ij, right), z_ij, which_i, k, n, mean, precision, ij) > y0 && jj++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  jj = 0;
  while(multinom_LL(replace_it(p_ij, ij, newx), z_ij, which_i, k, n, mean, precision, ij) < y0 && jj++ < ncontract) {
    if(newx < pij) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(jj == ncontract) return pij;
  return newx;
}

double cond_mv_mean_multinom(NumericVector x, IntegerVector refidx, NumericVector mean, NumericMatrix Q, int ij) {
  int i = refidx[ij];
  double mm = mean[i];
  for(int jj = 0; jj < x.size(); jj++) {
    int j = refidx[jj];
    if(ij != jj &&  j >= 0) {
      mm -= Q(i, j)*(x[jj] - mean[j]) / Q(i, i);
    }
  }
  return mm;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_multinom_mv(NumericMatrix p_ij, LogicalMatrix z_ij, IntegerVector which_i, LogicalVector is_ref, NumericMatrix k_ij, NumericMatrix n_ij, NumericMatrix mean, NumericMatrix Q, double w, int nexpand, int ncontract) {
  NumericMatrix out = clone(p_ij);
  IntegerVector refidx(is_ref.size());
  int ridx = 0;
  for(int ij = 0; ij < refidx.size(); ij++) {
    if(is_ref[ij]) {
      refidx[ij] = -1;
    } else {
      refidx[ij] = ridx++;
    }
  }

  for(int r=0; r < p_ij.nrow(); r++) {
    NumericVector kk = k_ij(r, _);
    NumericVector nn = n_ij(r, _);
    NumericVector mm = mean(r, _);
    for(int ij=0; ij < p_ij.ncol(); ij++) {
      if(is_ref[ij]) {
        out(r, ij) = 0.0;
        continue;
      }

      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean_multinom(out(r, _), refidx, mm, Q, ij);
      out(r, ij) = one_multinom_slice(out(r, _), z_ij(r, _), which_i, kk[ij], nn[ij], mmm, Q(refidx[ij], refidx[ij]), ij, w, nexpand, ncontract);
    }
  }
  return out;
}


