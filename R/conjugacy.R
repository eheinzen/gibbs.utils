
#' Draw Conjugate Samples
#'
#' Sample from the posterior (conditional on all other parameters) in the conjugate setting.
#'
#' @param y,x realizations from the distribution whose parameter is being drawn. For multivariate conjugacy,
#'   this is an n-by-p matrix
#' @param Q,tau the precision of the (multivariate) normal distribution from which \code{y} comes
#' @param mu0 the prior mean of \code{mu} or \code{beta}. The default (\code{NULL}) is the same as an appropriately-dimensioned
#'   vector/matrix of all zeros, but a little faster.
#' @param Q0,tau0 the prior precision of \code{mu}.
#' @param mult An optional multiplier for the prior precision; useful in some cases
#' @param mu the mean of the normal distribution from which \code{y} comes
#' @param a0,b0 the parameters (shape and rate) of the gamma distribution prior on \code{tau},
#'   or the parameters (shape1 and shape2) of the beta distribution prior on \code{p}.
#' @param a The shape parameter for the gamma distribution
#' @param k The number of successes
#' @param n The number of trials
#' @param V0,v0 the parameters (matrix and degrees of freedom) of the Wishart prior on \code{Q} or \code{V}
#' @param V,U the precision matrices for the matrix-normal distribution
#' @param X the data matrix on \code{beta}
#' @param beta the coefficients on \code{X}
#' @param XtX,XtU,XtUX,ytUy,V0_inv,Xbeta,newQ,newQ.chol pre-computed "shortcut" arguments for efficiency reasons
#' @param diag If \code{TRUE}, both \code{V} and \code{Q0} are assumed (but not confirmed!) to be diagonal,
#'   which can speed up the Cholesky decomposition up to 100x. Default \code{FALSE}.
#' @param zero A matrix of ones and zeros, the same size as the beta to sample. Zero indicates a structural zero in the beta.
#'   In the event that this is specified, everything returned is of size \code{sum(zero)}.
#' @param ... Other arguments. Examples include \code{verbose=}, \code{take.chol=}, and
#'   \code{Rstruct=}.
#' @param params.only Should just a list of the updated parameters be returned?
#' @name conjugacy
NULL

#' @rdname conjugacy
#' @export
conj_norm_mu <- function(y, tau, mu0 = 0, tau0 = 0.001, ..., mult = 1, params.only = FALSE) {
  n <- length(y)
  newtau <- mult*tau0 + n*tau
  mean <- (tau0*mu0 + tau*sum(y)) / newtau
  if(params.only) return(gu_params(mu = mean, sd = 1 / sqrt(newtau), tau = newtau))
  rnorm(
    1,
    mean = mean,
    sd = 1 / sqrt(newtau)
  )
}

#' @rdname conjugacy
#' @export
conj_mvnorm_mu <- function(y, Q, mu0 = NULL, Q0 = diag(0.001, p), ..., newQ = mult*Q0 + n*Q, newQ.chol = gu_chol(newQ), mult = 1, params.only = FALSE) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  p <- ncol(y)
  n <- nrow(y)
  Q0mu0 <- if(is.null(mu0)) 0 else Q0 %*% mu0
  b <- Q0mu0 + Q %*% colSums(y)

  if(params.only) {
    newQ.inv <- chol2inv(newQ.chol)
    mu <- drop(newQ.inv %*% b)
    gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv, b = drop(b))
  } else {
    newQ <- if(!missing(newQ.chol) && is.matrix(newQ.chol)) {
      gu_chol(newQ.chol, take.chol = FALSE) # the cholesky shouldn't be taken anymore after this
    } else {
      if(!missing(newQ.chol)) warning("'newQ.chol=' is being ignored")
      newQ
    }
    drop(spam::rmvnorm.canonical(1, b = b, Q = newQ, ...))
  }
}

#' @rdname conjugacy
#' @export
conj_matnorm_mu <- function(y, V, U = NULL, mu0 = NULL, Q0, ...,
                            newQ = V %x% U + Q0, newQ.chol = gu_chol(newQ), diag = FALSE, zero = NULL, params.only = FALSE) {
  if(!is.matrix(y)) stop("'y' must be a matrix")
  if(is.matrix(mu0)) mu0 <- vec(mu0)
  if(is.null(U)) U <- spam::diag.spam(1, nrow(Q0) / nrow(V))
  m <- ncol(U)
  p <- nrow(V)

  z <- if(is.null(zero)) m*p else {
    if(!is.matrix(zero) || nrow(zero) != m || ncol(zero) != p) stop("'zero' must be of dimension ", m, " x ", p)
    sum(zero == 1)
  }
  if(any(dim(Q0) != z)) stop("'Q0' should be of dimension ", z, " x ", z)
  if(!is.null(mu0) && length(mu0) != z) stop("'mu0' must be of length ", z)

  if(diag) {
    if(!missing(newQ) || !missing(newQ.chol)) warning("Arguments 'newQ' and 'newQ.chol' are being ignored because diag = TRUE")
    if(is.null(mu0)) mu0 <- 0
    Q0 <- spam::diag(Q0)
    if(!is.null(zero)) {
      Q0 <- vec(replace(zero, zero == 1, Q0))
      mu0 <- vec(replace(zero, zero == 1, mu0))
    }
    Q0.mat.lst <- asplit(matrix(Q0, nrow = m, ncol = p), 2)
    mu0.lst <- asplit(matrix(mu0, nrow = m, ncol = p), 2)
    Uy.lst <- asplit(U %*% y, 2)
    zero.lst <- if(is.null(zero)) rep_len(list(NULL), p) else asplit(zero == 1, 2)
    tmp <- Map(function(v, q, m, uyy, zz) {
      newQ <- v * U + spam::diag(q)
      b <- q*m + vec(uyy * v)
      if(!is.null(zz)) {
        newQ <- newQ[zz, zz, drop = FALSE]
        b <- b[zz]
      }

      if(params.only) {
        newQ.inv <- chol_inv(newQ)
        mu <- drop(newQ.inv %*% b)
        gu_params(mu = mu, Q.inv = newQ.inv, Q = newQ, b = drop(b))
      } else {
        drop(spam::rmvnorm.canonical(1, b = b, Q = newQ, ...))
      }
    }, v = spam::diag(V), q = Q0.mat.lst, m = mu0.lst, uyy = Uy.lst, zz = zero.lst)
    if(params.only) {
      return(gu_params(mu = lapply(tmp, "[[", "mu"), Q = lapply(tmp, "[[", "Q"), Q.inv = lapply(tmp, "[[", "Q.inv"), b = lapply(tmp, "[[", "b")))
    }
    do.call(c, tmp)

  } else {
    Q0mu0 <- if(is.null(mu0)) 0 else Q0 %*% mu0
    UyV <- U %*% y %*% V
    b <- Q0mu0 + if(is.null(zero)) vec(UyV) else UyV[zero == 1]
    if(missing(newQ) && missing(newQ.chol) && !is.null(zero)) {
      rz <- row(zero)[zero == 1]
      cz <- col(zero)[zero == 1]
      newQ <- matrix(0, z, z) #this method is faster than using rep, apparently
      ii <- as.vector(row(newQ))
      jj <- as.vector(col(newQ))
      newQ[] <- V[cbind(cz[ii], cz[jj])] * U[cbind(rz[ii], rz[jj])] + Q0
    }

    if(params.only) {
      newQ.inv <- chol2inv(newQ.chol)
      mu <- drop(newQ.inv %*% b)
      gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv, b = drop(b))
    } else {
      newQ <- if(!missing(newQ.chol) && is.matrix(newQ.chol)) {
        gu_chol(newQ.chol, take.chol = FALSE) # the cholesky shouldn't be taken anymore after this
      } else {
        if(!missing(newQ.chol)) warning("'newQ.chol=' is being ignored")
        newQ
      }
      drop(spam::rmvnorm.canonical(1, b = b, Q = newQ, ...))
    }
  }
}

#' @rdname conjugacy
#' @export
conj_lm_beta <- function(y, X, XtX = t(X) %*% X, tau, mu0 = NULL, Q0, ...,
                         newQ = tau * XtX + Q0, newQ.chol = gu_chol(newQ), params.only = FALSE) {
  Q0mu0 <- if(is.null(mu0)) 0 else Q0 %*% mu0
  b <- Q0mu0 + tau*t(X) %*% y

  if(params.only) {
    newQ.inv <- chol2inv(newQ.chol)
    mu <- drop(newQ.inv %*% b)
    gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv, b = drop(b))
  } else {
    newQ <- if(!missing(newQ.chol) && is.matrix(newQ.chol)) {
      gu_chol(newQ.chol, take.chol = FALSE) # the cholesky shouldn't be taken anymore after this
    } else {
      if(!missing(newQ.chol)) warning("'newQ.chol=' is being ignored")
      newQ
    }
    drop(spam::rmvnorm.canonical(1, b = b, Q = newQ, ...))
  }
}

#' @rdname conjugacy
#' @export
conj_matlm_beta <- function(y, X, V, U = NULL, mu0 = NULL, Q0, ..., XtU = if(is.null(U)) t(X) else t(X) %*% U, XtUX = XtU %*% X,
                            newQ = V %x% XtUX + Q0, newQ.chol = gu_chol(newQ), diag = FALSE, zero = NULL, params.only = FALSE) {
  if(!is.matrix(y)) stop("'y' must be a matrix")
  if(is.matrix(mu0)) mu0 <- vec(mu0)
  m <- ncol(XtUX)
  p <- nrow(V)

  z <- if(is.null(zero)) m*p else {
    if(!is.matrix(zero) || nrow(zero) != m || ncol(zero) != p) stop("'zero' must be of dimension ", m, " x ", p)
    sum(zero == 1)
  }
  if(any(dim(Q0) != z)) stop("'Q0' should be of dimension ", z, " x ", z)
  if(!is.null(mu0) && length(mu0) != z) stop("'mu0' must be of length ", z)

  if(diag) {
    if(!missing(newQ) || !missing(newQ.chol)) warning("Arguments 'newQ' and 'newQ.chol' are being ignored because diag = TRUE")
    if(is.null(mu0)) mu0 <- 0
    Q0 <- spam::diag(Q0)
    if(!is.null(zero)) {
      Q0 <- vec(replace(zero, zero == 1, Q0))
      mu0 <- vec(replace(zero, zero == 1, mu0))
    }
    Q0.mat.lst <- asplit(matrix(Q0, nrow = m, ncol = p), 2)
    mu0.lst <- asplit(matrix(mu0, nrow = m, ncol = p), 2)
    XtUy.lst <- asplit(XtU %*% y, 2)
    zero.lst <- if(is.null(zero)) rep_len(list(NULL), p) else asplit(zero == 1, 2)
    tmp <- Map(function(v, q, m, xtuyy, zz) {

      newQ <- v * XtUX + spam::diag(q)
      b <- q*m + vec(xtuyy * v)
      if(!is.null(zz)) {
        newQ <- newQ[zz, zz, drop = FALSE]
        b <- b[zz]
      }

      if(params.only) {
        newQ.inv <- chol_inv(newQ)
        mu <- drop(newQ.inv %*% b)
        gu_params(mu = mu, Q.inv = newQ.inv, Q = newQ, b = drop(b))
      } else {
        drop(spam::rmvnorm.canonical(1, b = b, Q = newQ, ...))
      }
    }, v = spam::diag(V), q = Q0.mat.lst, m = mu0.lst, xtuyy = XtUy.lst, zz = zero.lst)
    if(params.only) {
      return(gu_params(mu = lapply(tmp, "[[", "mu"), Q = lapply(tmp, "[[", "Q"), Q.inv = lapply(tmp, "[[", "Q.inv"), b = lapply(tmp, "[[", "b")))
    }
    do.call(c, tmp)

  } else {
    Q0mu0 <- if(is.null(mu0)) 0 else Q0 %*% mu0
    XtUyV <- XtU %*% y %*% V
    b <- Q0mu0 + if(is.null(zero)) vec(XtUyV) else XtUyV[zero == 1]
    if(missing(newQ) && missing(newQ.chol) && !is.null(zero)) {
      rz <- row(zero)[zero == 1]
      cz <- col(zero)[zero == 1]
      newQ <- matrix(0, z, z) #this method is faster than using rep, apparently
      ii <- as.vector(row(newQ))
      jj <- as.vector(col(newQ))
      newQ[] <- V[cbind(cz[ii], cz[jj])] * XtUX[cbind(rz[ii], rz[jj])] + Q0
    }

    if(params.only) {
      newQ.inv <- chol2inv(newQ.chol)
      mu <- drop(newQ.inv %*% b)
      gu_params(mu = mu, Q = newQ, Q.inv = newQ.inv, b = drop(b))
    } else {
      newQ <- if(!missing(newQ.chol) && is.matrix(newQ.chol)) {
        gu_chol(newQ.chol, take.chol = FALSE) # the cholesky shouldn't be taken anymore after this
      } else {
        if(!missing(newQ.chol)) warning("'newQ.chol=' is being ignored")
        newQ
      }
      drop(spam::rmvnorm.canonical(1, b = b, Q = newQ, ...))
    }
  }
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
conj_mvnorm_Q <- function(y, mu = NULL, V0, v0, V0_inv = chol_inv(V0), params.only = FALSE) {
  if(!is.matrix(y)) y <- matrix(y, nrow = 1)
  p <- ncol(y)
  n <- nrow(y)
  if(is.null(mu)) {
    ymu <- y
  } else if(!is.matrix(mu) && length(mu) %in% c(1, p)) {
    ymu <- y - matrix(mu, nrow = n, ncol = p, byrow = TRUE)
  } else if(is.matrix(mu)) {
    stopifnot(identical(dim(y), dim(mu)))
    ymu <- y - mu
  } else stop("'mu' must be a matrix of the same dimension as y or a vector with length 1 or ncol(y)")
  V2.inv <- V0_inv + t(ymu) %*% (ymu)
  V2 <- chol_inv(V2.inv)
  if(params.only) return(gu_params(V = V2, V.inv = V2.inv, v = n + v0))
  rWishart(
    1,
    Sigma = V2,
    df = n + v0
  )[, , 1]
}

#' @rdname conjugacy
#' @export
conj_matnorm_V <- function(y, mu = NULL, U = NULL, V0, v0, ..., ytUy = t(ymu) %*% U %*% ymu, V0_inv = chol_inv(V0), params.only = FALSE) {
  if(!is.matrix(y)) stop("'y' must be a matrix")
  if((missing(U) || is.null(U)) && missing(ytUy)) {
    return(conj_mvnorm_Q(y = y, mu = mu, v0 = v0, V0_inv = V0_inv, params.only = params.only))
  }
  n <- nrow(y)
  if(is.null(mu)) {
    ymu <- y
  } else {
    stopifnot(identical(dim(y), dim(mu)))
    ymu <- y - mu
  }
  V2.inv <- V0_inv + ytUy
  V2 <- chol_inv(V2.inv)
  if(params.only) return(gu_params(V = V2, V.inv = V2.inv, v = n + v0))
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



# non-normal --------------------------------------------------------------
#' @rdname conjugacy
#' @export
conj_binom_p <- function(k, n, a0 = 1, b0 = 1, params.only = FALSE) {
  a <- sum(k) + a0
  b <- sum(n) - sum(k) + b0
  if(params.only) return(gu_params(a = a, b = b))
  stats::rbeta(1, a, b)
}


#' @rdname conjugacy
#' @export
conj_gamma_b <- function(x, a, a0, b0, params.only = FALSE) {
  a <- a0 + length(x)
  b <- b0 + sum(x)
  if(params.only) return(gu_params(a = a, b = b))
  stats::rgamma(1, a, b)
}


