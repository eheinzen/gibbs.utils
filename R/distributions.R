
#' Statistical Distributions
#'
#' @param x,y vector of quantiles
#' @param mu the mean. The default (\code{NULL})- is the same as an appropriately-dimensioned
#'  vector/matrix of all zeros, but a little faster.
#' @param sd the standard deviation
#' @param tau,Q the precision (matrix)
#' @param byrow Should the densities be summed (FALSE, the default) or returned separately by row (TRUE).
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param detQ,detU,detV Pre-computed log-determinants
#' @param log Should the log-density be returned?
#' @param n Number of deviates to produce.
#' @param rate the Exponential rate parameter. If zero, a normal distribution is used instead. If negative,
#'   the problem is flipped and calculated using \code{-rate}.
#' @param a Bounds for the truncated exponential distribution.
#' @param b Bounds for the truncated exponential distribution or the multivariate normal canonical representation parameter.
#' @param q A quantile
#' @param p A probability
#' @param lower.tail Should lower tail probabilities be returned (default) or upper?
#' @param log.p Should the log be returned?
#' @param alpha the Dirichlet parameter vector or matrix.
#' @param meanlog,sdlog the parameters of the underlying log-normal distribution.
#' @details
#'   \code{dmvnorm_diff(x, y, mu, Q)} is equivalent to (but usually twice as fast as)
#'   \code{dmvnorm(x, mu, Q) - dmvnorm(y, mu, Q)}. Likewise \code{dmatnorm_diff}.
#' @name distributions
NULL

#' @rdname distributions
#' @export
dmvnorm <- function(x, mu = NULL, Q, detQ = determinant(Q, logarithm = TRUE)$modulus, log = TRUE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    if(!is.null(mu)) mu <- matrix(mu, nrow = 1)
  }
  p <- ncol(x)
  n <- nrow(x)

  xmu <- if(!is.null(mu)) {
    stopifnot(identical(dim(x), dim(mu)))
    x - mu
  } else x

  num <- n/2 * as.numeric(detQ)
  # tr(xmu Q xmu^T) = tr(xmu^T xmu Q) = vec(xmu)^T vec(xmu Q)
  num <- num - 0.5*sum(xmu * (xmu %*% Q))

  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}

#' @rdname distributions
#' @export
dmvlnorm <- function(x, mu = NULL, Q, detQ = determinant(Q, logarithm = TRUE)$modulus, log = TRUE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    if(!is.null(mu)) mu <- matrix(mu, nrow = 1)
  }
  p <- ncol(x)
  n <- nrow(x)

  lx <- log(x)
  xmu <- if(!is.null(mu)) {
    stopifnot(identical(dim(x), dim(mu)))
    lx - mu
  } else lx

  num <- n/2 * as.numeric(detQ)
  # tr(xmu Q xmu^T) = tr(xmu^T xmu Q) = vec(xmu)^T vec(xmu Q)
  num <- num + sum(-0.5*xmu * (xmu %*% Q) - lx)

  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}

#' @rdname distributions
#' @export
dmvnorm_diff <- function(x, y, mu, Q, b, log = TRUE, byrow = FALSE) {
  if((has.b <- missing(mu)) + missing(b) != 1) stop("Exactly one of 'mu' or 'b' must be specified")
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    y <- matrix(y, nrow = 1)
    if(has.b) {
      b <- matrix(b, nrow = 1)
    } else {
      mu <- matrix(mu, nrow = 1)
    }
  }

  FUN <- if(byrow) rowSums else sum
  if(has.b) {
    # tr((b^T x -  0.5*xQx^T) - (b^T y -  0.5*yQy^T))
    # tr(-0.5*(XQX^T - YQY^T)) + tr(b^T(x - y)))
    # tr(-0.5*(XQX^T - YQY^T)) + vec(b)^T vec(x - y)
    # -0.5*tr((X - Y)Q(X^T + Y^T)) + vec(b)^T vec(x - y)
    # -0.5*tr(Q(X^T + Y^T)(X - Y)) + vec(b)^T vec(x - y)
    # -0.5*vec((X + Y)Q)^T vec(X - Y) + vec(b)^T vec(x - y)
    # (-0.5*vec((X + Y)Q)^T + vec(b)^T)vec(X - Y)
    num <- FUN((-0.5*(x + y) %*% Q + b) * (x - y))
  } else {
    # tr((x - mu) Q (x - mu)^T - (y - mu) Q (y - mu)^T)
    # tr(XQX^T - YQY^T - 2(X - Y)Q mu^T)
    # tr((X - Y)Q(X^T + Y^T - 2 mu^T))
    # tr((X + Y - 2 mu)^T (X - Y) Q)
    num <- -0.5*FUN((x + y - 2*mu) * ((x - y) %*% Q))
  }
  if(log) num else exp(num)
}

#' @rdname distributions
#' @export
dmvlnorm_diff <- function(x, y, mu, Q, log = TRUE, byrow = FALSE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    y <- matrix(y, nrow = 1)
    mu <- matrix(mu, nrow = 1)
  }
  lx <- log(x)
  ly <- log(y)
  FUN <- if(byrow) rowSums else sum
  num <- FUN(-0.5*(lx + ly - 2*mu) * ((lx - ly) %*% Q) - (lx - ly))
  if(log) num else exp(num)
}

#' @rdname distributions
#' @export
dnorm_diff <- function(x, y, mu, tau, log = TRUE, byrow = FALSE) {
  FUN <- if(byrow) identity else sum
  num <- -0.5*FUN((x + y - 2*mu) * (x - y) * tau)
  if(log) num else exp(num)
}

#' @rdname distributions
#' @export
dlnorm_diff <- function(x, y, mu, tau, log = TRUE, byrow = FALSE) {
  FUN <- if(byrow) identity else sum
  lx <- log(x)
  ly <- log(y)
  num <- -0.5*FUN((lx + ly - 2*mu + 2/tau) * (lx - ly) * tau)
  if(log) num else exp(num)
}

#' @rdname distributions
#' @export
dmatnorm <- function(x, mu = NULL, V, U = NULL, detV = determinant(V, logarithm = TRUE)$modulus, detU = determinant(U, logarithm = TRUE)$modulus, log = TRUE) {
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(p == dim(V))
  Vdet <- as.numeric(detV)

  if(!is.null(U)) {
    stopifnot(n == dim(U))
    Udet <- as.numeric(detU)
  } else Udet <- 0

  xmu <- if(!is.null(mu)) {
    stopifnot(identical(dim(x), dim(mu)))
    x - mu
  } else x
  xV <- xmu %*% V
  UxV <- if(is.null(U)) xV else U %*% xV
  # tr(V x^T U x) = vec(x) vec(U x V)

  num <- 0.5*(p*Udet + n*Vdet - sum(xmu * UxV))
  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}

#' @rdname distributions
#' @export
dmatnorm_diff <- function(x, y, mu, V, U = NULL, log = TRUE) {
  xyV <- (x - y) %*% V
  UxV <- if(is.null(U)) xyV else U %*% xyV
  num <- -0.5*sum((x + y - 2*mu) * UxV)
  if(log) num else exp(num)
}



#' @rdname distributions
#' @export
dlogitnorm <- function(x, mu, sd, log = FALSE) {
  out <- stats::dnorm(logit(x), mean = mu, sd = sd, log = TRUE) - log(x) - log(1 - x)
  if(!log) exp(out) else out
}


#' @rdname distributions
#' @export
rdirich <- function(n, alpha) {
  K <- if(is.matrix(alpha)) ncol(alpha) else length(alpha)
  alpha <- matrix(alpha, nrow = n, ncol = K, byrow = TRUE)
  gam <- stats::rgamma(length(alpha), alpha, rate = 1)
  dim(gam) <- dim(alpha)
  gam / rowSums(gam)
}
