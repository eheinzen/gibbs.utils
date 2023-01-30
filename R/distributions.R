
#' Statistical Distributions
#'
#' @param x,y vector of quantiles
#' @param mu the mean
#' @param sd the standard deviation
#' @param tau,Q the precision (matrix)
#' @param byrow Should the densities be summed (FALSE, the default) or returned separately by row (TRUE).
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param detQ,detU,detV Pre-computed log-determinants
#' @param log Should the log-density be returned?
#' @param n Number of deviates to produce.
#' @param rate the Exponential rate parameter
#' @param a,b Bounds for the truncated exponential distribution.
#' @param q A quantile
#' @param p A probability
#' @param lower.tail Should lower tail probabilities be returned (default) or upper?
#' @param log.p Should the log be returned?
#' @param alpha the Dirichlet parameter vector or matrix.
#' @details
#'   \code{dmvnorm_diff(x, y, mu, Q)} is equivalent to (but usually twice as fast as)
#'   \code{dmvnorm(x, mu, Q) - dmvnorm(y, mu, Q)}. Likewise \code{dmatnorm_diff}.
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
dmvlnorm <- function(x, mu, Q, detQ = determinant(Q, logarithm = TRUE)$modulus, log = TRUE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    mu <- matrix(mu, nrow = 1)
  }
  stopifnot(identical(dim(x), dim(mu)))
  p <- ncol(x)
  n <- nrow(x)

  lx <- log(x)
  xmu <- lx - mu
  num <- n/2 * as.numeric(detQ)
  # tr(xmu Q xmu^T) = tr(xmu^T xmu Q) = vec(xmu)^T vec(xmu Q)
  num <- num + sum(-0.5*xmu * (xmu %*% Q) - lx)

  den <- n*p/2 * log(2*pi)
  if(log) num - den else exp(num - den)
}

#' @rdname distributions
#' @export
dmvnorm_diff <- function(x, y, mu, Q, log = TRUE, byrow = FALSE) {
  if(!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
    y <- matrix(y, nrow = 1)
    mu <- matrix(mu, nrow = 1)
  }
  # tr((x - mu) Q (x - mu)^T - (y - mu) Q (y - mu)^T)
  # tr(XQX^T - YQY^T - 2(X - Y)Q mu^T)
  # tr((X - Y)Q(X^T + Y^T - 2 mu^T))
  # tr((X + Y - 2 mu)^T (X - Y) Q)
  FUN <- if(byrow) rowSums else sum
  num <- -0.5*FUN((x + y - 2*mu) * ((x - y) %*% Q))
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

#' @rdname distributions
#' @export
rtruncexp <- function(n, rate = 1, a = 0, b = Inf) {
  stopifnot(b > a)
  u <- stats::runif(n)
  qtruncexp(p = u, rate = rate, a = a, b = b)
}

#' @rdname distributions
#' @export
dtruncexp <- function(x, rate = 1, a = 0, b = Inf, log = FALSE) {
  stopifnot(b > a)
  ll <- log(rate) + rate*a - rate*x - log(1 - exp(-rate*(b - a)))
  ll[x < a] <- -Inf
  ll[x > b] <- -Inf

  if(log) ll else exp(ll)
}

#' @rdname distributions
#' @export
ptruncexp <- function(q, rate = 1, a = 0, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(b > a)
  num <- exp(-rate*(q - a)) - 1
  den <- exp(-rate*(b - a)) - 1

  p <- num/den
  p[q > b] <- 1
  p[q < a] <- 0

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  p
}

#' @rdname distributions
#' @export
qtruncexp <- function(p, rate = 1, a = 0, b = Inf) {
  stopifnot(b > a, 0 <= p, p <= 1)
  a - log(p*exp(-rate*(b - a)) - p + 1)/rate
}
