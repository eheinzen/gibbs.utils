#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_pois_ratio(double L, double proposal, double k, double mean, double precision, int acceptance) {
  if(acceptance == 2) return 0.0;
  double ratio = pois_LL(proposal, k, mean, precision);
  ratio -= pois_LL(L, k, mean, precision);
  return ratio;
}


void qt_pois_approx(double around, double k, double mean, double tau, double& outmean, double& outsd) {
  // don't calculate this more than once
  double ep = exp(around);

  // -H(x)
  double newtau = tau + ep;
  outsd = 1.0/sqrt(newtau);

  // x - H^-1(x) g(x)
  double newmean = around + (k - ep - tau*(around - mean))/newtau;
  outmean = newmean;
}


void one_qt_pois_proposal_ratio(double L, double k, double mean, double precision, double& outproposal, double& outratio, int acceptance) {

  double prop_mean, prop_sd;
  qt_pois_approx(L, k, mean, precision, prop_mean, prop_sd);
  double proposal = R::rnorm(prop_mean, prop_sd);
  outproposal = proposal;

  double ratio = one_m_pois_ratio(L, proposal, k, mean, precision, acceptance);
  if(acceptance == 0) {
    double orig_mean, orig_sd;
    qt_pois_approx(proposal, k, mean, precision, orig_mean, orig_sd);
    ratio -= R::dnorm(proposal, prop_mean, prop_sd, 1);
    ratio += R::dnorm(L, orig_mean, orig_sd, 1);
  }
  outratio = ratio;
}

void gamma_pois_approx(double k, double mean, double tau, double& outalpha, double& outbeta) {
  double invtau = 1.0/tau;
  double m = exp(mean + 0.5*invtau);
  double mm = m*m;
  double v = (exp(invtau) - 1.0) * mm;
  double alpha = mm / v + k;
  double beta = m / v + 1.0;
  outalpha = alpha;
  outbeta = beta;
}

void one_gamma_pois_proposal_ratio(double L, double mult, double k, double mean, double precision, double& outproposal, double& outratio, int acceptance) {

  double alpha, beta;
  gamma_pois_approx(k, mean, precision, alpha, beta);
  double alpha2 = alpha/mult, scale2 = mult/beta;
  double proposal = R::rgamma(alpha2, scale2);
  double lproposal = log(proposal);
  outproposal = lproposal;

  double ratio = one_m_pois_ratio(L, lproposal, k, mean, precision, acceptance);
  ratio -= lproposal - L;
  if(acceptance == 0) {
    ratio -= (alpha2 - 1.0)*lproposal - proposal/scale2;
    ratio += (alpha2 - 1.0)*L - exp(L)/scale2;
  }
  outratio = ratio;
}



// [[Rcpp::export]]
NumericVector mh_pois(int method, NumericVector L, NumericVector proposal, NumericVector k, LogicalVector k_na,
                      NumericVector mean, NumericVector precision, int acceptance) {
  NumericVector out = clone(L);
  LogicalVector accept(L.size());

  for(int i=0; i < L.size(); i++) {
    if(k_na[i]) {
      out[i] = R::rnorm(mean[i], 1.0/sqrt(precision[i]));
      accept[i] = true;
      continue;
    }

    double ratio, prop;
    if(method == 1) { // 'proposal' is ignored
      one_qt_pois_proposal_ratio(L[i], k[i], mean[i], precision[i], prop, ratio, acceptance);
    } else if(method == 2) {
      one_gamma_pois_proposal_ratio(L[i], proposal[i], k[i], mean[i], precision[i], prop, ratio, acceptance);
    } else {
      prop = proposal[i];
      ratio = one_m_pois_ratio(L[i], prop, k[i], mean[i], precision[i], acceptance);
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
NumericVector mh_pois_mv(int method, NumericMatrix L, NumericMatrix proposal, NumericMatrix k, LogicalMatrix k_na,
                         NumericMatrix mean, NumericMatrix Q, LogicalVector use_norm, NumericMatrix norm, int acceptance) {
  NumericMatrix out = clone(L);
  LogicalMatrix accept(L.nrow(), L.ncol());

  int nrm = 0;
  for(int r=0; r < L.nrow(); r++) {
    if(use_norm[r]) {
      out(r, _) = norm(nrm++, _);
      for(int i=0; i < L.ncol(); i++) {
        accept(r, i) = true;
      }
      continue;
    }

    NumericVector kk = k(r, _);
    LogicalVector kk_na = k_na(r, _);
    NumericVector mm = mean(r, _);
    for(int i = 0; i < L.ncol(); i++) {
      // safe not to replace out(r, i) because it's not used in this calculation
      double mmm = cond_mv_mean(out(r, _), mm, Q, i);

      if(kk_na[i]) {
        out(r, i) = R::rnorm(mmm, 1.0/sqrt(Q(i, i)));
        accept(r, i) = true;
        continue;
      }

      double ratio, prop;
      if(method == 1) { // 'proposal' is ignored
        one_qt_pois_proposal_ratio(out(r, i), kk[i], mmm, Q(i, i), prop, ratio, acceptance);
      } else if(method == 2) {
        one_gamma_pois_proposal_ratio(out(r, i), proposal(r, i), kk[i], mmm, Q(i, i), prop, ratio, acceptance);
      } else {
        prop = proposal(r, i);
        ratio = one_m_pois_ratio(out(r, i), proposal(r, i), kk[i], mmm, Q(i, i), acceptance);
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
