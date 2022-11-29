#' Sample a multinomial regression rate
#'
#' @param p the previous iteration of the logit-probability, as an array of dimension \code{r x i x j},
#'   where \code{r} is the number of realizations, \code{i} is the number of (correlated)
#'   betas (that is, the dimension of \code{precision}), and \code{j} is the number of
#'   classes to split between. It is expected that \code{p[, , ref] == 0}. This can also be "flattened" to a
#'   matrix of dimension \code{r x (i x j)}.
#' @param z a matrix (or array) of zeros and ones. The zeros determine so-called "structural zeros": outcomes
#'   which are not possible. If two-dimensional, it is assumed not to change over the first dimension (\code{r}) of \code{p}.
#'   This is of dimension \code{i x j} or \code{r x i x j}.
#' @param zmax For advanced use: a matrix, usually computed from the last two dimensions of \code{z}.
#' @param k the realized value from the binomial distribution; the same size as \code{p}.
#' @param mean the prior mean for \code{p}. Must be either a scalar or an array of size \code{r x i x (j-1)} (when \code{p} is 3-D) or
#'   a "flattened" matrix of the same size.
#' @param precision an array of dimension \code{i x i x (j-1)} or a "flattened" matrix of the same size.
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
sample_multinom_reg <- function(p, z, k, mean, precision, method = c("slice"), ref = c("first", "last"), ...,
                                zmax = NULL, width = 1, nexpand = 10, ncontract = 100) {
  method <- match.arg(method)

  nr <- nrow(p)
  z <- TRUE & z
  if(length(dim(z)) == 3) {
    tmpzmax <- apply(z, 2:3, max)
    if(nrow(z) != nr) stop(glue::glue("The first dimension of 'z' ({nrow(z)}) must match the first dimension of 'p' ({nr})"))
    dim(z) <- c(nr, prod(dim(z)[2:3]))
  } else {
    tmpzmax <- z
    z <- matrix(z, nrow = nr, ncol = length(z), byrow = TRUE)
  }
  if(is.null(zmax)) {
    zmax <- tmpzmax
  } else if(!all(zmax >= tmpzmax)) {
    stop("'zmax' has fewer nonzero entries than 'z'")
  } else zmax <- TRUE & zmax

  if(is.character(ref)) {
    ref <- match.arg(ref)
    ref <- if(ref == "first") 1 else ncol(zmax)
  } else {
    if(!(ref %in% seq_len(ncol(zmax)))) stop(glue::glue("'ref' must be between 1 and {ncol(zmax))}"))
  }
  is_ref <- as.vector(col(zmax) == ref)
  which_i <- as.vector(row(zmax))

  d <- dim(p)
  if(length(d) == 3) {
    if(!identical(d[2:3], dim(zmax))) stop(glue::glue("The second and third dimensions of 'p' ({d[2]} x {d[3]}) must match 'z' ({nrow(zmax)} x {ncol(zmax)})"))
    if(!identical(d, dim(k))) stop(glue::glue("The dimensions of 'k' must match the dimensions of 'p' ({d[1]} x {d[2]} x {d[3]})"))

    dim(p) <- dim(k) <- c(nr, d[2]*d[3])
    dij <- d[2]*(d[3]-1L)

    dm <- dim(mean)
    if(identical(dm, c(nr, d[2], d[3]-1L))) {
      dim(mean) <- c(nr, dij)
    } else stop(glue::glue("The dimensions of 'mean' ({dm[1]} x {dm[2]} x {dm[3]}) should be {nr} x {d[2]} x {d[3]-1}"))

    dp <- dim(precision)
    if(identical(dp, c(d[2], d[2], d[3]-1L))) {
      prec_ij <- matrix(0, dij, dij)
      idx <- do.call(rbind, lapply(seq_len(d[3]-1), function(j) {
        cbind(
          rep(1:d[2], times = d[2]),
          rep(1:d[2], each = d[2])
        ) + (j-1)*d[2]
      }))
      prec_ij[idx] <- precision
      precision <- prec_ij
    } else stop(glue::glue("The dimensions of 'precision' ({dp[1]} x {dp[2]} x {dp[3]}) should be {d[2]} x {d[2]} x {d[3]-1}"))

  } else {
    if(d[2] == sum(zmax)) {
      zmax.vec <- as.vector(zmax)
      is_ref <- is_ref[zmax.vec]
      z <- z[, zmax.vec, drop = FALSE]
      which_i <- which_i[zmax.vec]
    } else {
      if(!identical(d[2], length(zmax))) stop(glue::glue("The second dimension of 'p' ({d[2]}) must match 'z' ({nrow(zmax)} x {ncol(zmax)})"))
    }
    if(!identical(d, dim(k))) stop(glue::glue("The dimensions of 'k' must match the dimensions of 'p' ({d[1]} x {d[2]})"))

    dm <- dim(mean)
    if(!identical(dm, c(nr, sum(!is_ref)))) {
      stop(glue::glue("The dimensions of 'mean' ({dm[1]} x {dm[2]}) should be {nr} x {sum(!is_ref)}"))
    }

    dp <- dim(precision)
    if(!identical(dp, c(sum(!is_ref), sum(!is_ref)))) {
      stop(glue::glue("The dimensions of 'precision' ({dp[1]} x {dp[2]}) should be {sum(!is_ref)} x {sum(!is_ref)}"))
    }
  }


  n <- vapply(seq_len(nrow(zmax)), function(i) rowSums(k[, which_i == i, drop = FALSE]), numeric(nr))
  dim(n) <- c(nr, nrow(zmax))
  n <- n[, which_i, drop = FALSE]
  stopifnot(
    identical(dim(p), dim(z)),
    identical(dim(p), dim(k)),
    identical(dim(p), dim(n)),
    length(which_i) == ncol(p),
    length(is_ref) == ncol(p),
    identical(dim(precision), c(sum(!is_ref), sum(!is_ref))),
    identical(dim(mean), c(nrow(p), sum(!is_ref)))

  )
  if(any(!z & k > 0)) stop("z == 0 but k > 0")
  if(any(p[, is_ref] != 0)) stop("p[, ref] != 0")

  p <- slice_sample_multinom_mv(
    p_ij = p,
    z_ij = z,
    k_ij = k,
    n_ij = n,
    which_i = which_i,
    is_ref = is_ref,
    mean = mean,
    Q = precision,
    w = width,
    nexpand = nexpand,
    ncontract = ncontract
  )
  dim(p) <- d
  p
}
