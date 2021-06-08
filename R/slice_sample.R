#' Slice sample a Poisson regression rate
#'
#' @param L the previous iteration of the log-rate
#' @param k the realized value from the Poisson distribution
#' @param mean the prior mean
#' @param precision the prior precision
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
ss_pois_reg <- function(L, k, mean, precision) {
  if(length(L) != length(k)) stop("'L' and 'k' must have the same length")
  mean <- check_one_or_all(mean, length(L))

  if(is.matrix(precision)) {
    slice_sample_pois_mv(L, k, mean, precision)
  } else {
    precision <- check_one_or_all(precision, length(L))
    slice_sample_pois(L, k, mean, precision)
  }
}


#' Slice sample a binomial regression rate
#'
#' @param p the previous iteration of the probability
#' @param k the realized value from the binomial distribution
#' @param n the number of trials
#' @param mean the prior mean
#' @param precision the prior precision
#' @details
#'   This function slice samples \code{p} conditional on \code{k}, \code{n}, \code{mean}, and \code{precision},
#'   where \code{k ~ Binom(n, p)} and \code{p ~ logitN(mean, precision)}.
#'
#'   This is vectorized over \code{p}, \code{k}, \code{n}, and \code{mean}.
#'
#'   The internals are defined in C++.
#' @seealso \code{\link{ss_pois_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling}
#' @export
ss_binom_reg <- function(p, k, mean, precision) {
  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  mean <- check_one_or_all(mean, length(p))

  precision <- check_one_or_all(precision, length(p))
  slice_sample_binom(L, k, n, mean, precision)
}
