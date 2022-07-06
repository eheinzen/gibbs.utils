
#' Statistical Distributions
#'
#' @param x vector of quantiles
#' @param mu the mean
#' @param sd the standard deviation
#' @param Q the precision matrix
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param detQ,detU,detV Pre-computed log-determinants
#' @param log Should the log-density be returned?
#' @name distributions
NULL

#' @rdname distributions
#' @export
dmvnorm <- function(x, mu, Q, detQ = determinant(Q, logarithm = TRUE)$modulus, log = TRUE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    mu <- matrix(mu, nrow = 1)
  }
  stopifnot(identical(dim(x), dim(mu)))
  p <- ncol(x)
  n <- nrow(x)

  xmu <- x - mu
  num <- n/2 * as.numeric(detQ)
  # tr(xmu Q xmu^T) = tr(xmu^T xmu Q) = vec(xmu)^T vec(xmu Q)
  num <- num - 0.5*sum(xmu * (xmu %*% Q))

  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}

#' @rdname distributions
#' @export
dmatnorm <- function(x, mu, V, U = NULL, detV = determinant(V, logarithm = TRUE)$modulus, detU = determinant(U, logarithm = TRUE)$modulus, log = TRUE) {
  stopifnot(identical(dim(x), dim(mu)))
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(p == dim(V))
  Vdet <- as.numeric(detV)

  if(!is.null(U)) {
    stopifnot(n == dim(U))
    Udet <- as.numeric(detU)
  } else Udet <- 0


  xmu <- x - mu
  xV <- xmu %*% V
  UxV <- if(is.null(U)) xV else U %*% xV
  # tr(V x^T U x) = vec(x) vec(U x V)

  num <- 0.5*(p*Udet + n*Vdet - sum(xmu * UxV))
  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}



#' @rdname distributions
#' @export
dlogitnorm <- function(x, mu, sd, log = FALSE) {
  out <- stats::dnorm(logit(x), mean = mu, sd = sd, log = TRUE) - log(x) - log(1 - x)
  if(!log) exp(out) else out
}

