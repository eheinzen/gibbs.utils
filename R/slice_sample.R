#' Slice sample a Poisson regression rate
#'
#' @param L the previous iteration of the log-rate
#' @param k the realized value from the Poisson distribution
#' @param mean the prior mean
#' @param precision the prior precision
#' @param ... Other arguments (not used)
#' @param w The width of each expansion (to the right and left each)
#' @param nexpand The maximum number of expansions (to the right and left each)
#' @param ncontract The maximum number of contractions. If this is exceeded, the original value is returned
#' @details
#'   This function slice samples \code{L} conditional on \code{k}, \code{mean}, and \code{precision},
#'   where \code{k ~ Pois(exp(L))} and \code{L ~ N(mean, precision)}.
#'
#'   This is vectorized over \code{L}, \code{k}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{L} is assumed to be multivariately distributed, and a different function is used.
#'
#'   The internals are defined in C++.
#' @seealso \code{\link{ss_binom_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling}
#' @export
ss_pois_reg <- function(L, k, mean, precision, ..., w = 1, nexpand = 10, ncontract = 100) {
  if(length(L) != length(k)) stop("'L' and 'k' must have the same length")
  mean <- check_one_or_all(mean, length(L))

  if(is.matrix(precision)) {
    FUN <- slice_sample_pois_mv
  } else {
    precision <- check_one_or_all(precision, length(L))
    FUN <- slice_sample_pois
  }
  FUN(L, k, mean, precision, w = w, nexpand = nexpand, ncontract = ncontract)
}


#' Slice sample a binomial regression rate
#'
#' @param p the previous iteration of the logit-probability
#' @param k the realized value from the binomial distribution
#' @param n the number of trials
#' @inheritParams ss_pois_reg
#' @details
#'   This function slice samples \code{p} conditional on \code{k}, \code{n}, \code{mean}, and \code{precision},
#'   where \code{k ~ Binom(n, expit(p))} and \code{p ~ N(mean, precision)}.
#'
#'   This is vectorized over \code{p}, \code{k}, \code{n}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{p} is assumed to be multivariately distributed, and a different function is used.
#'
#'   The internals are defined in C++.
#' @seealso \code{\link{ss_pois_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling}
#' @export
ss_binom_reg <- function(p, k, n, mean, precision, ..., w = 1, nexpand = 10, ncontract = 100) {
  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  mean <- check_one_or_all(mean, length(p))

  if(is.matrix(precision)) {
    FUN <- slice_sample_binom_mv
  } else {
    precision <- check_one_or_all(precision, length(p))
    FUN <- slice_sample_binom
  }
  FUN(p, k, n, mean, precision, w = w, nexpand = nexpand, ncontract = ncontract)
}
