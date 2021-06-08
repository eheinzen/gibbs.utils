
#' Statistical Distributions
#'
#' @param x vector of quantiles
#' @param mu the mean
#' @param sd the standard deviation
#' @param Q the precision matrix
#' @param V,U the precision matrices for the matrix-normal distribution, or (V) the parameter to the Wishart distribution
#' @param v the degree-of-freedom parameter for the Wishart distribution
#' @param log Should the log-density be returned?
#' @rdname distributions
#' @export
dmvnorm <- function(x, mu, Q, log = TRUE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    mu <- matrix(mu, nrow = 1)
  }
  stopifnot(identical(dim(x), dim(mu)))
  p <- ncol(x)
  n <- nrow(x)

  num <- n/2 * as.numeric(determinant(Q, logarithm = TRUE)$modulus)
  xmu <- x - mu
  for(i in seq_len(n)) {
    num <- num - 0.5*as.numeric(t(x[i, ] - mu[i, ]) %*% Q %*% (x[i, ] - mu[i, ]))
  }

  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}


#' @rdname distributions
#' @export
dmatnorm <- function(x, mu, V, U, log = TRUE) {
  stopifnot(identical(dim(x), dim(mu)))
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(p == dim(V), n == dim(U))

  Udet <- as.numeric(determinant(U, logarithm = TRUE)$modulus)
  Vdet <- as.numeric(determinant(V, logarithm = TRUE)$modulus)

  xmu <- as.numeric(x - mu)
  num <- 0.5*(p*Udet + n*Vdet) - 0.5*as.numeric(t(xmu) %*% (V %x% U) %*% xmu)
  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}


#' @rdname distributions
#' @export
dlogitnorm <- function(x, mu, sd, log = FALSE) {
  out <- stats::dnorm(logit(x), mean = mu, sd = sd, log = TRUE) - log(x) - log(1 - x)
  if(!log) exp(out) else out
}


#' @rdname distributions
#' @export
rwish <- function(V, v){
  if (v < nrow(V)) stop("v is less than the dimension of V in rwish().\n")

  p <- nrow(V)
  L <- chol(V)
  A <- matrix(0, p, p)
  diag(A) <- sqrt(stats::rchisq(p, v:(v - p + 1)))
  if(p > 1) {
    A[upper.tri(A)] <- stats::rnorm(p*(p-1)/2)
  }
  crossprod(A %*% L)
}


