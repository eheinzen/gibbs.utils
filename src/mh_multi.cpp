#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_multinom_ratio(const double p, const double proposal,
                            const double sum_exp_p, const bool z_ij,
                            const double k, const double n,
                            const double mean, const double precision,
                            const int acceptance) {
  if(acceptance == 2) return 0.0;
  return multinom_LL(proposal, sum_exp_p, z_ij, k, n, mean, precision) -
    multinom_LL(p, sum_exp_p, z_ij, k, n, mean, precision);
}


void qt_multinom_approx(const double around, const double sum_exp_p,
                        const double k, const double n,
                        const double mean, const double tau,
                        double& outmean, double& outsd) {
  // don't calculate these more than once
  double ep = exp(around);
  double ep1 = sum_exp_p + ep;
  double ep2 = ep1*ep1;

  // -H(x)
  double newtau = tau + n*ep*sum_exp_p / ep2;
  outsd = 1.0/sqrt(newtau);

  // x - H^-1(x) g(x)
  double newmean = around + (k - n*ep/ep1 - tau*(around - mean))/newtau;
  outmean = newmean;
}

void one_qt_multinom_proposal_ratio(const double p, const double sum_exp_p, const bool z_ij,
                                    const double k, const double n,
                                    const double mean, const double precision,
                                    double& outproposal, double& outratio,
                                    const int acceptance) {

  double prop_mean, prop_sd;
  qt_multinom_approx(p, sum_exp_p, /*z_ij,*/ k, n, mean, precision, prop_mean, prop_sd);
  double proposal = R::rnorm(prop_mean, prop_sd);
  outproposal = proposal;

  double ratio = one_m_multinom_ratio(p, proposal, sum_exp_p, z_ij, k, n, mean, precision, acceptance);
  if(acceptance == 0) {
    double orig_mean, orig_sd;
    qt_multinom_approx(proposal, sum_exp_p, /*z_ij,*/ k, n, mean, precision, orig_mean, orig_sd);
    ratio -= R::dnorm(proposal, prop_mean, prop_sd, 1);
    ratio += R::dnorm(p, orig_mean, orig_sd, 1);
  }
  outratio = ratio;
}


// [[Rcpp::export]]
NumericVector mh_multinom_mv(const bool qt, const NumericMatrix p_ij,const  NumericMatrix proposal,
                             const LogicalMatrix z_ij, const IntegerVector which_i, const LogicalVector is_ref,
                             const NumericMatrix k_ij, const NumericMatrix n_ij,
                             const NumericMatrix mean, const NumericMatrix Q,
                             const LogicalVector use_norm, const NumericMatrix norm,
                             const bool diag, const int acceptance) {
  NumericMatrix out = clone(p_ij);
  LogicalMatrix accept(p_ij.nrow(), p_ij.ncol());
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
      for(int ij=0; ij < p_ij.ncol(); ij++) {
        accept(r, ij) = true;
      }
      continue;
    }
    NumericVector kk = k_ij(r, _);
    NumericVector nn = n_ij(r, _);
    NumericVector mm = mean(r, _);
    LogicalVector zz = z_ij(r, _);
    NumericVector tmp = out(r, _);
    NumericVector propp = proposal(r, _);
    LogicalVector acc = accept(r, _);

    for(int ij=0; ij < p_ij.ncol(); ij++) {
      if(is_ref[ij]) {
        tmp[ij] = 0.0;
        acc[ij] = true;
        continue;
      }
      // safe not to replace tmp[ij] because it's not used in this calculation
      double mmm;
      if(diag) {
        mmm = mm[refidx[ij]];
      } else {
        mmm = cond_mv_mean_multinom(tmp, refidx, mm, Q, ij);
      }

      if(nn[ij] == 0.0) {
        tmp[ij] = R::rnorm(mmm, 1.0/sqrt(Q(refidx[ij], refidx[ij])));
        acc[ij] = true;
        continue;
      }

      int i = which_i[ij];
      double sum_exp_p = 0.0;
      for(int iijj = 0; iijj < zz.size(); iijj++) {
        if(iijj != ij && zz[iijj] && which_i[iijj] == i) {
          sum_exp_p += exp(tmp[iijj]);
        }
      }

      double ratio, prop;
      if(qt) { // 'proposal' is ignored
        one_qt_multinom_proposal_ratio(tmp[ij], sum_exp_p, zz[ij], kk[ij], nn[ij], mmm, Q(refidx[ij], refidx[ij]), prop, ratio, acceptance);
      } else {
        prop = propp[ij];
        ratio = one_m_multinom_ratio(tmp[ij], prop, sum_exp_p, zz[ij], kk[ij], nn[ij], mmm, Q(refidx[ij], refidx[ij]), acceptance);
      }
      bool a = accept_reject(ratio);
      acc[ij] = a;
      if(a) {
        tmp[ij] = prop;
      }
    }
    out(r, _) = tmp;
    accept(r, _) = acc;
  }
  out.attr("accept") = accept;
  return out;
}
