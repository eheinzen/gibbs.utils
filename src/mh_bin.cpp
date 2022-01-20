#include "utils.h"
#include <Rcpp.h>
using namespace Rcpp;

double one_m_ratio(double p, double proposal, double k, double n, double mean, double precision) {
  double ratio = binom_LL(proposal, k, n, mean, precision);
  ratio -= binom_LL(p, k, n, mean, precision);
  return ratio;
}

bool accept_reject(double ratio) {
  return log(R::runif(0.0, 1.0)) <= ratio;
}

// [[Rcpp::export]]
NumericVector m_binom(NumericVector p, NumericVector proposal, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision,
                      LogicalVector use_norm, NumericVector norm) {
  NumericVector out = clone(p);
  LogicalVector accept(p.size());

  int nrm = 0;
  for(int i=0; i < p.size(); i++) {
    if(use_norm[i]) {
      out[i]= norm[nrm++];
      accept[i] = true;
      continue;
    }

    double ratio = one_m_ratio(p[i], proposal[i], k[i], n[i], mean[i], precision[i]);
    bool a = accept_reject(ratio);
    accept[i] = a;
    if(a) {
      out[i] = proposal[i];
    }
  }
  out.attr("accept") = accept;
  return out;
}



// [[Rcpp::export]]
NumericVector m_binom_mv(NumericMatrix p, NumericMatrix proposal, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q,
                         LogicalVector use_norm, NumericMatrix norm) {
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
        out(r, i) = R::rnorm(mmm, 1/sqrt(Q(i, i)));
        accept(r, i) = true;
        continue;
      }
      double ratio = one_m_ratio(out(r, i), proposal(r, i), kk[i], nn[i], mmm, Q(i, i));
      bool a = accept_reject(ratio);
      accept(r, i) = a;
      if(a) {
        out(r, i) = proposal(r, i);
      }
    }
  }
  out.attr("accept") = accept;
  return out;
}



// [[Rcpp::export]]
NumericVector qt_binom_approx(double around, double k, double n, double mean, double tau) {
  // don't calculate these more than once
  double ep = exp(around);
  double ep1 = 1.0 + ep;
  double ep2 = ep1*ep1;

  // -H(x)
  double newtau = tau + n*ep / ep2;

  // x - H^-1(x) g(x)
  double newmean = around + (k + tau*mean - tau*around - n*ep/ep1)/newtau;

  return NumericVector::create(newmean, 1.0/sqrt(newtau));
}


NumericVector one_qt_proposal_ratio(double p, double k, double n, double mean, double precision) {

  NumericVector prop_params = qt_binom_approx(p, k, n, mean, precision);
  double proposal = R::rnorm(prop_params[0], prop_params[1]);

  NumericVector orig_params = qt_binom_approx(proposal, k, n, mean, precision);

  double ratio = binom_LL(proposal, k, n, mean, precision);
  ratio -= R::dnorm(proposal, prop_params[0], prop_params[1], 1);

  ratio -= binom_LL(p, k, n, mean, precision);
  ratio += R::dnorm(p, orig_params[0], orig_params[1], 1);

  return NumericVector::create(proposal, ratio);
}


// [[Rcpp::export]]
NumericVector qt_binom(NumericVector p, NumericVector k, NumericVector n, NumericVector mean, NumericVector precision,
                       LogicalVector use_norm, NumericVector norm) {
  NumericVector out = clone(p);
  LogicalVector accept(p.size());

  int nrm = 0;
  for(int i=0; i < p.size(); i++) {
    if(use_norm[i]) {
      out[i]= norm[nrm++];
      accept[i] = true;
      continue;
    }

    NumericVector prop_ratio = one_qt_proposal_ratio(p[i], k[i], n[i], mean[i], precision[i]);
    bool a = accept_reject(prop_ratio[1]);
    accept[i] = a;
    if(a) {
      out[i] = prop_ratio[0];
    }
  }
  out.attr("accept") = accept;
  return out;
}

// [[Rcpp::export]]
NumericVector qt_binom_mv(NumericMatrix p, NumericMatrix k, NumericMatrix n, NumericMatrix mean, NumericMatrix Q,
                          LogicalVector use_norm, NumericMatrix norm) {
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
        out(r, i) = R::rnorm(mmm, 1/sqrt(Q(i, i)));
        accept(r, i) = true;
        continue;
      }
      NumericVector prop_ratio = one_qt_proposal_ratio(out(r, i), kk[i], nn[i], mmm, Q(i, i));
      bool a = accept_reject(prop_ratio[1]);
      accept(r, i) = a;
      if(a) {
        out(r, i) = prop_ratio[0];
      }
    }
  }
  out.attr("accept") = accept;
  return out;
}
