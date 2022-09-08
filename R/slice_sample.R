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
#'   In the case that \code{k} is \code{NA} (akin to \code{\link{ss_binom_reg}} in the case \code{n == k == 0}),
#'   slice sampling is ignored in favor of a normal draw.
#'   In the special case when an entire (multivariate) row of \code{k} is \code{NA}, the entire row
#'   is simultaneously drawn (in R); when only some elements are \code{NA}, they're drawn univariately (in C++).
#'   Note that only the latter is used when \code{proposal} indicates multivariate MH.
#'
#'   This is vectorized over \code{L}, \code{k}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{L} is assumed to be multivariately distributed, and a different function is used.
#'
#'   The internals are defined in C++.
#' @seealso \code{\link{ss_binom_reg}}, \code{\link{ss_multinom_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling},
#' \url{https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm}, \url{https://arxiv.org/pdf/1308.0657.pdf}
#' @name pois_reg
NULL

#' @rdname pois_reg
#' @export
ss_pois_reg <- function(L, k, mean, precision, ..., w = 1, nexpand = 10, ncontract = 100) {
  if(length(L) != length(k)) stop("'L' and 'k' must have the same length")
  mean <- check_one_or_all(mean, length(L))
  d <- dim(L)

  if(is.matrix(precision)) {
    tmp <- is.matrix(L) + is.matrix(k)
    if(!(tmp %in% c(0, 2))) stop("Both of 'L' and 'k' must be matrices, or neither must be.")
    if(!is.matrix(L)) {
      L <- matrix(L, nrow = 1)
      k <- matrix(k, nrow = 1)
    }
    dim(mean) <- dim(L)
    use_norm <- rowSums(!is.na(k)) == 0
    norm <- if(any(use_norm)) {
      mean[use_norm, , drop = FALSE] + chol_mvrnorm(sum(use_norm), mu = 0, Precision = precision)
    } else matrix()
    out <- slice_sample_pois_mv(L, k, is.na(k), mean, precision, use_norm = use_norm, norm = norm, w = w, nexpand = nexpand, ncontract = ncontract)
  } else {
    precision <- check_one_or_all(precision, length(L))
    out <- slice_sample_pois(L, k, is.na(k), mean, precision, w = w, nexpand = nexpand, ncontract = ncontract)
  }

  dim(out) <- d # could be NULL
  out
}


#' Slice sample a binomial regression rate
#'
#' @param p the previous iteration of the logit-probability
#' @param k the realized value from the binomial distribution
#' @param n the number of trials
#' @inheritParams ss_pois_reg
#' @inherit ss_pois_reg seealso
#' @details
#'   \code{ss_binom_reg} slice samples and \code{mh_binom_reg} Metropolis-samples
#'   \code{p} conditional on \code{k}, \code{n}, \code{mean}, and \code{precision},
#'   where \code{k ~ Binom(n, expit(p))} and \code{p ~ N(mean, precision)}.
#'
#'   In the case that \code{n} is zero, slice sampling is ignored in favor of a normal draw.
#'   In the special case when an entire (multivariate) row of \code{n} is zero, the entire row
#'   is simultaneously drawn (in R); when only some elements are zero, they're drawn univariately (in C++).
#'   Note that only the latter is used when \code{proposal} indicates multivariate MH.
#'
#'   Both vectorized over \code{p}, \code{k}, \code{n}, and \code{mean}. If \code{precision} is a matrix,
#'   \code{p} is assumed to be multivariately distributed, and different internals are used.
#'
#'   The internals are defined in C++.
#' @name binom_reg
NULL

#' @rdname binom_reg
#' @export
ss_binom_reg <- function(p, k, n, mean, precision, ..., w = 1, nexpand = 10, ncontract = 100) {
  if(length(p) != length(k) || length(p) != length(n)) stop("'p' and 'k' and 'n' must all have the same length")
  stopifnot(k <= n)
  mean <- check_one_or_all(mean, length(p))
  d <- dim(p)

  if(is.matrix(precision)) {
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
    out <- slice_sample_binom_mv(p, k, n, mean, precision, use_norm = use_norm, norm = norm, w = w, nexpand = nexpand, ncontract = ncontract)
  } else {
    precision <- check_one_or_all(precision, length(p))
    out <- slice_sample_binom(p, k, n, mean, precision, w = w, nexpand = nexpand, ncontract = ncontract)
  }
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
#' @param ref One of \code{"first"} or \code{"last"} or a numeric value, denoting which "column" (the third dimension)
#'    of \code{p} is the reference.
#'   If \code{"first"}, then \code{mean} and \code{precision} map to the second through j-th elements of \code{p}.
#'   If \code{"last"}, then \code{mean} and \code{precision} map to the first through (j-1)-th
#'   elements of \code{p}.
#'   If \code{2}, then \code{mean} and \code{precision} map to the first and third through j-th elements of \code{p}.
#'   Etc.
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

  if(is.character(ref)) {
    ref <- match.arg(ref)
    ref <- if(ref == "first") 1 else d[3]
  } else {
    if(!(ref %in% seq_len(d[3]))) stop("'ref' must be between 1 and the length of the third dimension of 'p'")
  }
  if(any(p[, , ref] != 0)) stop("p[, , ref] != 0")

  if(length(mean) == 1) mean <- array(mean, dim = d2)
  if(!identical(dim(mean), d2)) {
    stop("'mean' must either be a scalar or have dim = r x i x (j-1): ", d2[1], " x ", d2[2], " x ", d2[3])
  }

  dp <- dim(precision)
  if(!identical(dp, d2[c(2, 2, 3)])) {
    stop("'precision' must be have dim = i x i x (j-1): ", d2[2], " x ", d2[2], " x ", d2[3])
  }

  err <- paste0("'z' must be a matrix with dim = i x j (", d[2], " x ", d[3], ") or an array with dim = r x i x j (", d[1], " x ", d[2], " x ", d[3], ")")
  if(length(dim(z)) == 3) {
    if(!identical(dim(z), d)) stop(err)
    z <- TRUE & z
    if(any(!z & k > 0)) stop("z == 0 but k > 0")
    z <- asplit(z, 1)
  } else {
    if(!identical(dim(z), d[2:3])) stop(err)
    z <- TRUE & z
    if(any(!z & colSums(k) > 0)) stop("z == 0 but k > 0")
  }
  n <- rowSums(k, dims = 2)
  subset_third <- function(x, i) {
    # this is in case 'x' has other dims of length 1
    out <- x[, , i]
    dim(out) <- dim(x)[1:2]
    out
  }
  FUN <- if(is.list(z)) slice_sample_multinom_mv_zlist else slice_sample_multinom_mv
  for(j in seq_len(d[3])[-ref]) {
    p[, , j] <- FUN(
      p_j = unclass(asplit(p, 1)),
      z = z,
      k = subset_third(k, j),
      n = n,
      p_i = subset_third(p, j),
      mean = subset_third(mean, j - (j > ref)),
      Q = subset_third(precision, j - (j > ref)),
      j = j - 1L,
      w = w,
      nexpand = nexpand,
      ncontract = ncontract
    )
  }
  p
}
