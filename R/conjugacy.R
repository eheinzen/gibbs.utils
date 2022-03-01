
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
#' @param a0,b0 the parameters (shape and rate) of the gamma distribution prior on \code{tau},
#'   or the parameters (shape1 and shape2) of the beta distribution prior on \code{p}.
#' @param k The number of successes
#' @param n The number of trials
#' @param V0,v0 the parameters (matrix and degrees of freedom) of the Wishart prior on \code{Q} or \code{V}
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param X the data matrix on \code{beta}
#' @param beta the coefficients on \code{X}
#' @param XtX,V0_inv,Xbeta,newQ.inv,A pre-computed "shortcut" arguments for efficiency reasons
#' @param diag If \code{TRUE}, both \code{V} and \code{Q0} are assumed (but not confirmed!) to be diagonal,
#'   which can speed up the Cholesky decomposition up to 100x. Default \code{FALSE}.
#' @param use.chol Logical, determining whether to use \code{MASS::\link{mvrnorm}()} (\code{FALSE}, the default) or \code{\link{chol_mvrnorm}}.
#' @param params.only Should just a list of the updated parameters be returned?
#' @rdname conjugacy
#' @export
conj_norm_mu <- function(y, tau, mu0 = 0, tau0 = 0.001, mult = 1, params.only = FALSE) {
  n <- length(y)
  newtau <- mult*tau0 + n*tau
  mean <- (tau0*mu0 + tau*sum(y)) / newtau
  if(params.only) return(gu_params(mu = mean, sd = 1 / sqrt(newtau)))
  rnorm(
    1,
    mean = mean,
    sd = 1 / sqrt(newtau)
  )
}

#' @rdname conjugacy
#' @export
conj_mvnorm_mu <- function(y, Q, mu0 = rep_len(0, p), Q0 = diag(0.001, p), newQ.inv = chol_inv(mult*Q0 + n*Q),
                           A = t(chol(newQ.inv)), mult = 1, use.chol = FALSE, params.only = FALSE) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  p <- ncol(y)
  n <- nrow(y)
  mu <- newQ.inv %*% (Q0 %*% mu0 + Q %*% colSums(y))
  if(params.only) return(gu_params(mu = mu, Q.inv = newQ.inv))
  if(use.chol) chol_mvrnorm(1, mu = mu, A = A) else MASS::mvrnorm(1, mu = mu, Sigma = newQ.inv)
}

#' @rdname conjugacy
#' @export
conj_matnorm_mu <- function(y, V, U = NULL, mu0, Q0, newQ.inv = chol_inv(V %x% U + Q0),
                            A = t(chol(newQ.inv)), diag = FALSE, use.chol = FALSE, params.only = FALSE) {
  if(!is.matrix(y)) stop("'y' must be a matrix")
  if(is.matrix(mu0)) mu0 <- as.numeric(mu0)
  if(is.null(U)) U <- diag(nrow(Q0) / nrow(V))

  if(diag) {
    if(!missing(newQ.inv) || !missing(A)) warning("Arguments 'newQ.inv' and 'A' are being ignored because diag = TRUE")
    Q0.mat.lst <- asplit(matrix(diag(Q0), nrow = ncol(X), ncol = nrow(V)), 2)
    mu0.lst <- asplit(matrix(mu0, nrow = ncol(X), ncol = nrow(V)), 2)
    y.lst <- asplit(y, 2)
    tmp <- Map(function(v, q, m, yy) {
      newQ <- v * U + diag(q)
      newQ.inv <- chol_inv(newQ)
      mu <- drop(newQ.inv %*% (q*m + as.numeric(U %*% yy * v)))
      if(params.only) return(gu_params(mu = mu, Q.inv = newQ.inv))
      FUN <- if(use.chol) chol_mvrnorm else MASS::mvrnorm
      FUN(1, mu = mu, Sigma = newQ.inv)
    }, v = diag(V), q = Q0.mat.lst, m = mu0.lst, yy = y.lst)
    if(params.only) {
      return(gu_params(mu = lapply(tmp, "[[", "mu"), Q = lapply(tmp, "[[", "Q"), Q.inv = lapply(tmp, "[[", "Q.inv")))
    }
    do.call(c, tmp)

  } else {
    mu <- drop(newQ.inv %*% (Q0 %*% mu0 + as.numeric(U %*% y %*% V)))
    if(params.only) return(gu_params(mu = mu, Q.inv = newQ.inv))
    if(use.chol) chol_mvrnorm(1, mu = mu, A = A) else MASS::mvrnorm(1, mu = mu, Sigma = newQ.inv)
  }
}

#' @rdname conjugacy
#' @export
conj_lm_beta <- function(y, X, XtX = crossprod(X), tau, mu0, Q0, newQ.inv = chol_inv(tau * XtX + Q0),
                         A = t(chol(newQ.inv)), use.chol = FALSE, params.only = FALSE) {
  mu <- newQ.inv %*% (Q0 %*% mu0 + tau*t(X) %*% y)
  if(params.only) return(gu_params(mu = mu, Q.inv = newQ.inv))
  if(use.chol) chol_mvrnorm(1, mu = mu, A = A) else MASS::mvrnorm(1, mu = mu, Sigma = newQ.inv)
}

#' @rdname conjugacy
#' @export
conj_matlm_beta <- function(y, X, V, U = NULL, mu0, Q0, diag = FALSE, use.chol = FALSE, params.only = FALSE) {
  if(!is.matrix(y)) stop("'y' must be a matrix")
  if(is.matrix(mu0)) mu0 <- as.numeric(mu0)

  XtU <- if(is.null(U)) t(X) else crossprod(X, U)
  XtUX <- XtU %*% X
  FUN <- if(use.chol) chol_mvrnorm else MASS::mvrnorm

  if(diag) {
    Q0.mat.lst <- asplit(matrix(diag(Q0), nrow = ncol(X), ncol = nrow(V)), 2)
    mu0.lst <- asplit(matrix(mu0, nrow = ncol(X), ncol = nrow(V)), 2)
    y.lst <- asplit(y, 2)
    tmp <- Map(function(v, q, m, yy) {
      newQ <- v * XtUX + diag(q)
      newQ.inv <- chol_inv(newQ)
      mu <- drop(newQ.inv %*% (q*m + as.numeric(XtU %*% yy * v)))
      if(params.only) return(gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv))
      FUN(1, mu = mu, Sigma = newQ.inv)
    }, v = diag(V), q = Q0.mat.lst, m = mu0.lst, yy = y.lst)
    if(params.only) {
      return(gu_params(mu = lapply(tmp, "[[", "mu"), Q = lapply(tmp, "[[", "Q"), Q.inv = lapply(tmp, "[[", "Q.inv")))
    }
    do.call(c, tmp)

  } else {
    newQ <- V %x% XtUX + Q0
    newQ.inv <- chol_inv(newQ)
    mu <- drop(newQ.inv %*% (Q0 %*% mu0 + as.numeric(XtU %*% y %*% V)))
    if(params.only) return(gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv))
    FUN(1, mu = mu, Sigma = newQ.inv)
  }

}

#' @rdname conjugacy
#' @export
conj_diagmatlm_beta <- function(y, X, V, U = NULL, mu0, Q0, use.chol = FALSE, params.only = FALSE) {
  if(!is.matrix(y)) stop("'y' must be a matrix")
  if(is.matrix(mu0)) mu0 <- as.numeric(mu0)

  p <- ncol(y)
  E <- replace(matrix(0, p^2, p), vapply(1:p, function(i) (i-1)*p^2 + (i-1)*p + i, NA_real_), 1)
  XtU <- if(is.null(U)) t(X) else crossprod(X, U)
  newQ <- crossprod(E, V %x% (XtU %*% X)) %*% E + Q0
  newQ.inv <- chol_inv(newQ)
  mu <- newQ.inv %*% (Q0 %*% mu0 + crossprod(E, as.numeric(XtU %*% y %*% V)))
  if(params.only) return(gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv))
  FUN <- if(use.chol) chol_mvrnorm else MASS::mvrnorm
  FUN(1, mu = mu, Sigma = newQ.inv)
}







# precisions --------------------------------------------------------------

#' @rdname conjugacy
#' @export
conj_norm_tau <- function(y, mu, a0 = 0.001, b0 = 0.001, params.only = FALSE) {
  a <- a0 + 0.5*length(y)
  b <- b0 + 0.5*sum((y - mu)^2)
  if(params.only) return(gu_params(a = a, b = b))
  rgamma(1, a, b)
}

#' @rdname conjugacy
#' @export
conj_mvnorm_Q <- function(y, mu, V0, v0, V0_inv = chol_inv(V0), params.only = FALSE) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  p <- ncol(y)
  n <- nrow(y)
  if(!is.matrix(mu) && length(mu) == p) {
    mu <- matrix(mu, nrow = n, ncol = p, byrow = TRUE)
  }
  stopifnot(identical(dim(y), dim(mu)))
  tmp <- crossprod(y - mu)
  V2 <- chol_inv(V0_inv + tmp)
  if(params.only) return(gu_params(V = V2, v = n + v0))
  rWishart(
    1,
    Sigma = V2,
    df = n + v0
  )[, , 1]
}

#' @rdname conjugacy
#' @export
conj_matnorm_V <- function(y, mu, U = NULL, V0, v0, V0_inv = chol_inv(V0), params.only = FALSE) {
  if(!is.matrix(y) || !is.matrix(mu)) stop("'y' and 'mu' must be matrices")
  n <- nrow(mu)
  stopifnot(identical(dim(y), dim(mu)))
  ymu <- y - mu
  ymu2 <- if(is.null(U)) t(ymu) else crossprod(ymu, U)

  V2 <- chol_inv(V0_inv + ymu2 %*% ymu)
  if(params.only) return(gu_params(V = V2, v = n + v0))
  rWishart(
    1,
    Sigma = V2,
    df = n + v0
  )[, , 1]
}

#' @rdname conjugacy
#' @export
conj_lm_tau <- function(y, X, beta, Xbeta = X %*% beta, a0 = 0.001, b0 = 0.001, params.only = FALSE) {
  conj_norm_tau(y = y, mu = Xbeta, a0 = a0, b0 = b0, params.only = params.only)
}

#' @rdname conjugacy
#' @export
conj_matlm_sigma <- function(y, X, beta, Xbeta = X %*% beta, U = NULL, V0, v0, V0_inv = chol_inv(V0), params.only = FALSE) {
  conj_matnorm_V(y = y, mu = Xbeta, U = U, v0 = v0, V0_inv = V0_inv, params.only = params.only)
}





# non-normal --------------------------------------------------------------
#' @rdname conjugacy
#' @export
conj_binom_p <- function(k, n, a0 = 1, b0 = 1, params.only = FALSE) {
  a <- sum(k) + a0
  b <- sum(n) - sum(k) + b0
  if(params.only) return(gu_params(a = a, b = b))
  stats::rbeta(1, a, b)
}





