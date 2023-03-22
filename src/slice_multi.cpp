#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_multinom_slice(const double p_ij, const double sum_exp_p, const bool z_ij,
                          const double k, const double n,
                          const double mean, const double precision,
                          const double w, const int nexpand, const int ncontract) {
  if(n == 0.0 || !z_ij) {
    return R::rnorm(mean, 1/sqrt(precision));
  }
  double y0 = multinom_LL(p_ij, sum_exp_p, z_ij, k, n, mean, precision) - R::rexp(1.0);
  double left = p_ij - w;
  double right = p_ij + w;
  int jj = 0;
  while(multinom_LL(left, sum_exp_p, z_ij, k, n, mean, precision) > y0 && jj++ < nexpand) {
    left -= w;
  }
  jj = 0;
  while(multinom_LL(right, sum_exp_p, z_ij, k, n, mean, precision) > y0 && jj++ < nexpand) {
    right += w;
  }
  double newx = R::runif(left, right);
  jj = 0;
  while(multinom_LL(newx, sum_exp_p, z_ij, k, n, mean, precision) < y0 && jj++ < ncontract) {
    if(newx < p_ij) {
      left = newx;
    } else {
      right = newx;
    }
    newx = R::runif(left, right);
  }
  if(jj == ncontract) return p_ij;
  return newx;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_multinom_mv(const NumericMatrix p_ij, const LogicalMatrix z_ij,
                                       const IntegerVector which_i, const LogicalVector is_ref,
                                       const NumericMatrix k_ij, const NumericMatrix n_ij,
                                       const NumericMatrix mean, const NumericMatrix Q,
                                       const bool diag,
                                       const LogicalVector use_norm, const NumericMatrix norm,
                                       const double w, const int nexpand, const int ncontract) {
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
  int nrm = 0;
  for(int r=0; r < p_ij.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      continue;
    }
    NumericVector kk = k_ij(r, _);
    NumericVector nn = n_ij(r, _);
    NumericVector mm = mean(r, _);
    LogicalVector zz = z_ij(r, _);
    NumericVector tmp = out(r, _);
    for(int ij=0; ij < p_ij.ncol(); ij++) {
      if(is_ref[ij]) {
        tmp[ij] = 0.0;
        continue;
      }

      // safe not to replace tmp[ij] because it's not used in this calculation
      double mmm;
      if(diag) {
        mmm = mm[refidx[ij]];
      } else {
        mmm = cond_mv_mean_multinom(tmp, refidx, mm, Q, ij);
      }

      int i = which_i[ij];
      double sum_exp_p = 0.0;
      for(int iijj = 0; iijj < zz.size(); iijj++) {
        if(iijj != ij && zz[iijj] && which_i[iijj] == i) {
          sum_exp_p += exp(tmp[iijj]);
        }
      }
      tmp[ij] = one_multinom_slice(tmp[ij], sum_exp_p, zz[ij], kk[ij], nn[ij], mmm, Q(refidx[ij], refidx[ij]), w, nexpand, ncontract);
    }
    out(r, _) = tmp;
  }
  return out;
}
