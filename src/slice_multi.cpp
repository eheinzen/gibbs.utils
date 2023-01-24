#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double one_multinom_slice(double p_ij, double sum_exp_p, bool z_ij, double k, double n, double mean, double precision, double w, int nexpand, int ncontract) {
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

double cond_mv_mean_multinom(NumericVector x, IntegerVector refidx, NumericVector mean, NumericMatrix Q, int ij) {
  int i = refidx[ij];
  double mm = mean[i];
  double Qii = Q(i, i);
  for(int jj = 0; jj < x.size(); jj++) {
    int j = refidx[jj];
    if(ij != jj &&  j >= 0) {
      mm -= Q(i, j)*(x[jj] - mean[j]) / Qii;
    }
  }
  return mm;
}

// [[Rcpp::export]]
NumericMatrix slice_sample_multinom_mv(NumericMatrix p_ij, LogicalMatrix z_ij, IntegerVector which_i, LogicalVector is_ref, NumericMatrix k_ij, NumericMatrix n_ij, NumericMatrix mean, NumericMatrix Q, bool diag, double w, int nexpand, int ncontract) {
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
    LogicalVector zz = z_ij(r, _);
    for(int ij=0; ij < p_ij.ncol(); ij++) {
      if(is_ref[ij]) {
        out(r, ij) = 0.0;
        continue;
      }

      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm;
      if(diag) {
        mmm = mm[refidx[ij]];
      } else {
        mmm = cond_mv_mean_multinom(out(r, _), refidx, mm, Q, ij);
      }

      int i = which_i[ij];
      double sum_exp_p = 0.0;
      for(int iijj = 0; iijj < zz.size(); iijj++) {
        if(iijj != ij && zz[iijj] && which_i[iijj] == i) {
          sum_exp_p += exp(out(r, iijj));
        }
      }
      out(r, ij) = one_multinom_slice(out(r, ij), sum_exp_p, zz[ij], kk[ij], nn[ij], mmm, Q(refidx[ij], refidx[ij]), w, nexpand, ncontract);
    }
  }
  return out;
}


