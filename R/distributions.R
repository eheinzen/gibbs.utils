
#' Statistical Distributions
#'
#' @param x vector of quantiles
#' @param mu the mean
#' @param sd the standard deviation
#' @param Q the precision matrix
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param log Should the log-density be returned?
#' @param use_trace Should the matrix trace be used, or a for-loop?
#' @name distributions
NULL

#' @rdname distributions
#' @export
dmvnorm <- function(x, mu, Q, log = TRUE, use_trace = FALSE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    mu <- matrix(mu, nrow = 1)
  }
  stopifnot(identical(dim(x), dim(mu)))
  p <- ncol(x)
  n <- nrow(x)

  xmu <- x - mu
  num <- n/2 * as.numeric(determinant(Q, logarithm = TRUE)$modulus)
  if(use_trace) {
    tmp <- xmu %*% tcrossprod(Q, xmu)
    num <- num - 0.5*sum(diag(tmp))
  } else {
    for(i in seq_len(n)) {
      num <- num - 0.5*as.numeric(xmu[i, , drop = FALSE] %*% tcrossprod(Q, xmu[i, , drop = FALSE]))
    }
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

