#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_binom_ratio(double p, double proposal, double k, double n, double mean, double precision, int acceptance) {
  if(acceptance == 2) return 0.0;
  return binom_LL(proposal, k, n, mean, precision) - binom_LL(p, k, n, mean, precision);
}


void qt_binom_approx(double around, double k, double n, double mean, double tau, double& outmean, double& outsd) {
  // don't calculate these more than once
  double ep = exp(around);
  double ep1 = 1.0 + ep;
  double ep2 = ep1*ep1;

  // -H(x)
  double newtau = tau + n*ep / ep2;
  outsd = 1.0/sqrt(newtau);

  // x - H^-1(x) g(x)
  double newmean = around + (k - n*ep/ep1 - tau*(around - mean))/newtau;
  outmean = newmean;
}

void one_qt_binom_proposal_ratio(double p, double k, double n, double mean, double precision, double& outproposal, double& outratio, int acceptance) {

  double prop_mean, prop_sd;
  qt_binom_approx(p, k, n, mean, precision, prop_mean, prop_sd);
  double proposal = R::rnorm(prop_mean, prop_sd);
  outproposal = proposal;

  double ratio = one_m_binom_ratio(p, proposal, k, n, mean, precision, acceptance);
  if(acceptance == 0) {
    double orig_mean, orig_sd;
    qt_binom_approx(proposal, k, n, mean, precision, orig_mean, orig_sd);
    ratio -= R::dnorm(proposal, prop_mean, prop_sd, 1);
    ratio += R::dnorm(p, orig_mean, orig_sd, 1);
  }
  outratio = ratio;
}


// [[Rcpp::export]]
NumericVector mh_binom(bool qt, NumericVector p, NumericVector proposal, NumericVector k, NumericVector n,
                       NumericVector mean, NumericVector precision, int acceptance) {
  NumericVector out = clone(p);
  LogicalVector accept(p.size());

  for(int i=0; i < p.size(); i++) {
    if(n[i] == 0.0) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      accept[i] = true;
      continue;
    }

    double ratio, prop;
    if(qt) { // 'proposal' is ignored
      one_qt_binom_proposal_ratio(p[i], k[i], n[i], mean[i], precision[i], prop, ratio, acceptance);
    } else {
      prop = proposal[i];
      ratio = one_m_binom_ratio(p[i], prop, k[i], n[i], mean[i], precision[i], acceptance);
    }
    bool a = accept_reject(ratio);
    accept[i] = a;
    if(a) {
      out[i] = prop;
    }
  }
  out.attr("accept") = accept;
  return out;
}


// [[Rcpp::export]]
NumericVector mh_binom_mv(bool qt, NumericMatrix p, NumericMatrix proposal, NumericMatrix k, NumericMatrix n,
                          NumericMatrix mean, NumericMatrix Q, LogicalVector use_norm, NumericMatrix norm, int acceptance) {
  NumericMatrix out = clone(p);
  LogicalMatrix accept(p.nrow(), p.ncol());

  int nrm = 0;
  for(int r=0; r < p.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      for(int i=0; i < p.ncol(); i++) {
        accept(r, i) = true;
      }
      continue;
    }

    NumericVector kk = k(r, _);
    NumericVector nn = n(r, _);
    NumericVector mm = mean(r, _);
    for(int i = 0; i < p.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);

      if(nn[i] == 0.0) {
        out(r, i) = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        accept(r, i) = true;
        continue;
      }

      double ratio, prop;
      if(qt) { // 'proposal' is ignored
        one_qt_binom_proposal_ratio(out(r, i), kk[i], nn[i], mmm, Q(i, i), prop, ratio, acceptance);
      } else {
        prop = proposal(r, i);
        ratio = one_m_binom_ratio(out(r, i), prop, kk[i], nn[i], mmm, Q(i, i), acceptance);
      }
      bool a = accept_reject(ratio);
      accept(r, i) = a;
      if(a) {
        out(r, i) = prop;
      }
    }
  }
  out.attr("accept") = accept;
  return out;
}
