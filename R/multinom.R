#' Sample a multinomial regression rate
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
#' @inheritParams sample_pois_reg
#' @details
#'   The internals are defined in C++.
#'
#'   In the case that \code{n} is zero or \code{z} is zero, slice sampling is ignored in favor of a normal draw.
#' @seealso \code{\link{sample_pois_reg}}, \code{\link{sample_binom_reg}}, \url{https://en.wikipedia.org/wiki/Slice_sampling}
#' @export
sample_multinom_reg <- function(p, z, k, mean, precision, method = c("slice"), ref = c("first", "last"), ..., width = 1, nexpand = 10, ncontract = 100) {
  method <- match.arg(method)

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
      w = width,
      nexpand = nexpand,
      ncontract = ncontract
    )
  }
  p
}
