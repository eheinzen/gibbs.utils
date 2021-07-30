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
  d <- dim(L)

  if(is.matrix(precision)) {
    FUN <- slice_sample_pois_mv
    tmp <- is.matrix(L) + is.matrix(k)
    if(!(tmp %in% c(0, 2))) stop("Both of 'L' and 'k' must be matrices, or neither must be.")
    if(!is.matrix(L)) {
      L <- matrix(L, nrow = 1)
      k <- matrix(k, nrow = 1)
    }
    dim(mean) <- dim(L)
  } else {
    precision <- check_one_or_all(precision, length(L))
    FUN <- slice_sample_pois
  }
  out <- FUN(L, k, mean, precision, w = w, nexpand = nexpand, ncontract = ncontract)
  dim(out) <- d # could be NULL
  out
}


#' Slice sample a binomial regression rate
#'
#' @param p the previous iteration of the logit-probability
#' @param k the realized value from the binomial distribution
#' @param n the number of trials
#' @inheritParams ss_pois_reg
#' @details
#'   \code{ss_binom_reg} slice samples and \code{mh_binom_reg} Metropolis-samples
#'   \code{p} conditional on \code{k}, \code{n}, \code{mean}, and \code{precision},
#'   where \code{k ~ Binom(n, expit(p))} and \code{p ~ N(mean, precision)}.
#'
#'   Both vectorized over \code{p}, \code{k}, \code{n}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{p} is assumed to be multivariately distributed, and different internals are used.
#'
#'   The internals are defined in C++.
#' @seealso \code{\link{ss_pois_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling},
#' \url{https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm}
#' @name binom_reg
#' @export
ss_binom_reg <- function(p, k, n, mean, precision, ..., w = 1, nexpand = 10, ncontract = 100) {
  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  mean <- check_one_or_all(mean, length(p))
  d <- dim(p)

  if(is.matrix(precision)) {
    FUN <- slice_sample_binom_mv
    tmp <- is.matrix(p) + is.matrix(k) + is.matrix(n)
    if(!(tmp %in% c(0, 3))) stop("All of 'p', 'k', and 'n' must be matrices, or none must be.")
    if(!is.matrix(p)) {
      p <- matrix(p, nrow = 1)
      k <- matrix(k, nrow = 1)
      n <- matrix(n, nrow = 1)
    }
    dim(mean) <- dim(p)
  } else {
    precision <- check_one_or_all(precision, length(p))
    FUN <- slice_sample_binom
  }
  out <- FUN(p, k, n, mean, precision, w = w, nexpand = nexpand, ncontract = ncontract)
  dim(out) <- d # could be NULL
  out
}
