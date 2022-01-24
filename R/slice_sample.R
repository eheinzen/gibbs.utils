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
#' @seealso \code{\link{ss_binom_reg}}, \code{\link{ss_multinom_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling}
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
#'   In the case that \code{n} is zero, slice sampling is ignored in favor of a normal draw.
#'   In the special case when an entire (multivariate) row of \code{n} is zero, the entire row
#'   is simultaneously drawn; when only some elements are zero, they're drawn univariately.
#'
#'   Both vectorized over \code{p}, \code{k}, \code{n}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{p} is assumed to be multivariately distributed, and different internals are used.
#'
#'   The internals are defined in C++.
#' @seealso \code{\link{ss_pois_reg}}, \code{\link{ss_multinom_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling},
#' \url{https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm}
#' @name binom_reg
#' @export
ss_binom_reg <- function(p, k, n, mean, precision, ..., w = 1, nexpand = 10, ncontract = 100) {
  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  stopifnot(k <= n)
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
    use_norm <- rowSums(n) == 0
    norm <- if(any(use_norm)) {
      mean[use_norm, , drop = FALSE] + chol_mvrnorm(sum(use_norm), mu = 0, Precision = precision)
    } else matrix()

  } else {
    precision <- check_one_or_all(precision, length(p))
    FUN <- slice_sample_binom
    use_norm <- n == 0
    norm <- if(any(use_norm)) {
      mean[use_norm] + stats::rnorm(sum(use_norm), mean = 0, sd = 1 / sqrt(precision[use_norm]))
    } else numeric()
  }
  out <- FUN(p, k, n, mean, precision, use_norm = use_norm, norm = norm, w = w, nexpand = nexpand, ncontract = ncontract)
  dim(out) <- d # could be NULL
  out
}



#' Slice sample a multinomial regression rate
#'
#' @param p the previous iteration of the logit-probability, as an array of dimension \code{r x i x j},
#'   where \code{r} is the number of realizations, \code{i} is the number of (correlated)
#'   betas (that is, the dimension of \code{precision}), and \code{j} is the number of
#'   classes to split between. It is expected that \code{p[, , ref] == 0}
#' @param z a matrix of zeros and ones. The zeros determine so-called "structural zeros": outcomes
#'   which are not possible. This is assumed not to change over the first dimension (\code{r}) of \code{p}.
#'   This is of dimension \code{i x j}.
#' @param k the realized value from the binomial distribution; the same size as \code{p} (\code{r x i x j}).
#' @param mean the prior mean for \code{p}. Must be either a scalar or an array of size \code{r x i x (j-1)}
#' @param precision an array of dimension \code{i x i x (j-1)}.
#'   The \code{j}-dimension is assumed to be independent.
#' @param ref One of \code{"first"} or \code{"last"}, denoting which column of \code{p} is the reference.
#'   If \code{"first"}, then \code{mean} and \code{precision} map to the second through j-th elements of \code{p}.
#'   If \code{"last"}, then \code{mean} and \code{precision} map to the first through (j-1)-th
#'   elements of \code{p}.
#' @inheritParams ss_pois_reg
#' @details
#'   The internals are defined in C++.
#'
#'   In the case that \code{n} is zero or \code{z} is zero, slice sampling is ignored in favor of a normal draw.
#' @seealso \code{\link{ss_pois_reg}}, \code{\link{ss_binom_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling}
#' @name multinom_reg
#' @export
ss_multinom_reg <- function(p, z, k, mean, precision, ref = c("first", "last"), ..., w = 1, nexpand = 10, ncontract = 100) {
  d <- dim(p)
  d2 <- replace(d, 3, d[3]-1L)
  if(length(d) != 3 || !identical(d, dim(k))) stop("'p' and 'k' must have the same (3D) dimensions")

  ref <- match.arg(ref)
  ref <- if(ref == "first") 1 else d[3]
  if(any(p[, , ref] != 0)) stop("p[, , ref] != 0")

  if(length(mean) == 1) mean <- array(mean, dim = d2)
  if(!identical(dim(mean), d2)) {
    stop("'mean' must either be a scalar or have dim = r x i x (j-1): ", d2[1], " x ", d2[2], " x ", d2[3])
  }

  dp <- dim(precision)
  if(!identical(dp, d2[c(2, 2, 3)])) {
    stop("'precision' must be have dim = i x i x (j-1): ", d2[2], " x ", d2[2], " x ", d2[3])
  }

  if(!identical(dim(z), d[2:3])) stop("'z' must be a matrix with dim = i x j")
  z <- TRUE & z
  if(any(!z & colSums(k))) stop("z == 0 but k > 0")
  n <- rowSums(k, dims = 2)
  subset_third <- function(x, i) {
    out <- x[, , i]
    dim(out) <- dim(x)[1:2]
    out
  }
  for(j in seq_len(d[3])[-ref]) {
    p[, , j] <- slice_sample_multinom_mv(
      p_j = unclass(asplit(p, 1)),
      z = z,
      k = subset_third(k, j),
      n = n,
      p_i = subset_third(p, j),
      mean = subset_third(mean, j - (ref == 1)),
      Q = subset_third(precision, j - (ref == 1)),
      j = j - 1L,
      w = w,
      nexpand = nexpand,
      ncontract = ncontract
    )
  }
  p
}
