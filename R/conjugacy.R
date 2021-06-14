
#' Draw Conjugate Samples
#'
#' Sample from the posterior (conditional on all other parameters) in the conjugate setting.
#'
#' @param y realizations from the distribution whose parameter is being drawn. For multivariate conjugacy,
#'   this is an n-by-p matrix
#' @param Q,tau the precision of the (multivariate) normal distribution from which \code{y} comes
#' @param mu0 the prior mean of \code{mu} or \code{beta}.
#' @param Q0,tau0 the prior precision of \code{mu}.
#' @param mult An optional multiplier for the prior precision; useful in some cases
#' @param mu the mean of the normal distribution from which \code{y} comes
#' @param a0,b0 the parameters (shape and rate) of the gamma distribution prior on \code{tau}.
#' @param V0,v0 the parameters (matrix and degrees of freedom) of the Wishart prior on \code{Q} or \code{V}
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param X the data matrix on \code{beta}
#' @param beta the coefficients on \code{X}
#' @param XtX,V0_inv,Xbeta pre-computed "shortcut" arguments for efficiency reasons
#' @rdname conjugacy
#' @export
conj_norm_mu <- function(y, tau, mu0 = 0, tau0 = 0.001, mult = 1) {
  n <- length(y)
  newtau <- mult*tau0 + n*tau
  stats::rnorm(
    1,
    mean = (tau0*mu0 + tau*sum(y)) / newtau,
    sd = 1 / sqrt(newtau)
  )
}

#' @rdname conjugacy
#' @export
conj_mvnorm_mu <- function(y, Q, mu0 = rep_len(0, p), Q0 = diag(0.001, p), mult = 1) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  p <- ncol(y)
  n <- nrow(y)
  newQ <- mult*Q0 + n*Q
  newQ.inv <- chol_inv(newQ)
  MASS::mvrnorm(
    1,
    mu = newQ.inv %*% (Q0 %*% mu0 + Q %*% colSums(y)),
    Sigma = newQ.inv
  )
}

#' @rdname conjugacy
#' @export
conj_matnorm_mu <- function(y, V, U = NULL, mu0, Q0) {
  if(is.matrix(y)) y <- as.numeric(y)
  if(is.null(U)) U <- diag(nrow(Q0) / nrow(V))
  VxU <- V %x% U
  newQ <- VxU + Q0
  newQ.inv <- chol_inv(newQ)
  MASS::mvrnorm(
    1,
    mu = newQ.inv %*% (Q0 %*% mu0 + VxU %*% y),
    Sigma = newQ.inv
  )
}

#' @rdname conjugacy
#' @export
conj_lm_beta <- function(y, X, XtX = crossprod(X), tau, mu0, Q0) {
  newQ <- tau * XtX + Q0
  newQ.inv <- chol_inv(newQ)
  MASS::mvrnorm(
    1,
    mu = newQ.inv %*% (Q0 %*% mu0 + tau*t(X) %*% y),
    Sigma = newQ.inv
  )
}

#' @rdname conjugacy
#' @export
conj_matlm_beta <- function(y, X, V, U = NULL, mu0, Q0) {
  if(is.matrix(y)) y <- as.numeric(y)
  if(is.matrix(mu0)) mu0 <- as.numeric(mu0)

  XtU <- if(is.null(U)) t(X) else crossprod(X, U)
  newQ <- V %x% (XtU %*% X) + Q0
  newQ.inv <- chol_inv(newQ)
  MASS::mvrnorm(
    1,
    mu = newQ.inv %*% (Q0 %*% mu0 + (V %x% XtU) %*% y),
    Sigma = newQ.inv
  )
}













#' @rdname conjugacy
#' @export
conj_norm_tau <- function(y, mu, a0 = 0.001, b0 = 0.001) {
  stats::rgamma(
    1,
    a0 + 0.5*length(y),
    b0 + 0.5*sum((y - mu)^2)
  )
}

#' @rdname conjugacy
#' @export
conj_mvnorm_Q <- function(y, mu, V0, v0, V0_inv = chol_inv(V0)) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  p <- ncol(y)
  n <- nrow(y)
  if(!is.matrix(mu) && length(mu) == p) {
    mu <- matrix(rep(mu, each = n), nrow = n, ncol = p)
  }
  stopifnot(identical(dim(y), dim(mu)))
  xmu <- y - mu
  tmp <- Reduce(`+`, lapply(1:n, function(i) {
    txmu <- xmu[i, , drop = FALSE]
    t(txmu) %*% txmu
  }))
  V2 <- chol_inv(V0_inv + tmp)
  rwish(
    V = V2,
    v = n + v0
  )
}

#' @rdname conjugacy
#' @export
conj_matnorm_V <- function(y, mu, U = NULL, V0, v0, V0_inv = chol_inv(V0)) {
  if(!is.matrix(y) || !is.matrix(mu)) stop("'y' and 'mu' must be matrices")
  n <- nrow(mu)
  stopifnot(identical(dim(y), dim(mu)))
  ymu <- y - mu
  ymu2 <- if(is.null(U)) t(ymu) else crossprod(ymu, U)

  V2 <- chol_inv(V0_inv + ymu2 %*% ymu)
  rwish(
    V = V2,
    v = n + v0
  )
}

#' @rdname conjugacy
#' @export
conj_lm_tau <- function(y, X, beta, Xbeta = X %*% beta, a0 = 0.001, b0 = 0.001) {
  conj_norm_tau(y = y, mu = Xbeta, a0 = a0, b0 = b0)
}

#' @rdname conjugacy
#' @export
conj_matlm_sigma <- function(y, X, beta, Xbeta = X %*% beta, U = NULL, V0, v0, V0_inv = chol_inv(V0)) {
  conj_matnorm_V(y = y, mu = Xbeta, U = U, v0 = v0, V0_inv = V0_inv)
}








