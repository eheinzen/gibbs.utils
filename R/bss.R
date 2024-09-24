
#' Calculate Bayesian Smoothing Spline (BSS) for Integers
#'
#' @param x for \code{bss}, an integer vector to calculate the basis expansion for.
#'   For \code{B1}, \code{B2}, \code{B4}: a numeric vector between 0 and 1.
#' @param n_basis Number of basis terms to use
#' @param first_two Should the first two BSS terms (linear and quadratic) be included?
#' @param from,to The min and max values that \code{x} can attain.
#' @references Curtis B. Storlie et al. (2015) Calibration of Computational Models
#'  With Categorical Parameters and Correlated Outputs via Bayesian Smoothing
#'  Spline ANOVA, Journal of the American Statistical Association, 110:509, 68-82,
#'  DOI: 10.1080/01621459.2014.979993
#' @details
#' Currently this does no linear interpolation and so relies on integer (discrete)
#' inputs.
#' @name bss
NULL


#' @rdname bss
#' @export
B1 <- function(x) x - 0.5

#' @rdname bss
#' @export
B2 <- function(x) x^2 - x + 1/6

#' @rdname bss
#' @export
B4 <- function(x) x^4 - 2*x^3 + x^2 - 1/30

#' @rdname bss
#' @export
bss <- function(x, n_basis, first_two = TRUE, from = min(x), to = max(x)) {
  if(!is.integer(x)) stop("This function only supports integer `x` at this point")
  if(anyNA(x)) stop("NAs in `x` are not supported")
  if(n_basis <= 0) stop("`n_basis` must be > 0")
  if(first_two && n_basis < 2) stop("`n_basis` must be >= 2 if `first_two` is `TRUE`")
  u <- from:to
  if(!all(x %in% u)) stop("Not all of `x` are in `from:to`.")
  mult <- if(length(u) < 360) ceiling(360/length(u)) else 1
  N <- mult*length(u)
  s <- seq(from = 0, to = 1, length.out = N+1)[-(N+1)]
  n_basis <- n_basis - 2*first_two
  if(n_basis > 0) {
    K365 <- matrix(-1/24*B4(abs(rep(s, times = N) - rep(s, each = N))), N, N)
    e365 <- eigen(K365, symmetric = TRUE)
    x365 <- e365$vectors[, 1:n_basis, drop = FALSE] %*% diag(sqrt(e365$values[1:n_basis]))
  } else x365 <- NULL
  if(first_two) {
    x365 <- cbind(B1(s), B2(s), x365)
  }

  idx <- mult*(x - from) + 1
  stopifnot(idx >= 1, idx <= nrow(x365))
  x365[idx, ]
}

