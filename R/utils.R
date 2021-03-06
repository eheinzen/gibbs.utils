
#' @useDynLib gibbs.utils
#' @importFrom Rcpp sourceCpp
NULL

check_one_or_all <- function(x, len) {
  if(length(x) == 1) rep_len(x, len) else if(length(x) == len) x else stop("'x' must be of length 1 or ", len)
}

#' Other Utilities
#'
#' @param x A matrix for \code{chol_inv} or a numeric vector for \code{expit}
#' @param p A numeric vector
#' @param n Number of samples to draw
#' @param mu,Sigma,Precision parameters for a multivariate normal distribution. One of \code{Sigma} and
#'   \code{Precision} must be supplied.
#' @param A the "square root" by which we multiply the independent normals. Used only for speed.
#' @param ... Other arguments (not in use at this time)
#' @param verbose Logical, indicating whether to print some stuff.
#' @param take.chol Should the cholesky be taken?
#' @details
#'   Unlike \code{MASS::\link[MASS]{mvrnorm}()}, \code{chol_mvrnorm()} uses the Cholesky decomposition
#'   to draw random samples. This is usually faster; the largest gains are seen when \code{n} is small
#'   (when \code{n} is large, the matrix multiplication becomes the limiting factor).
#' @name utilities
NULL


#' @rdname utilities
#' @export
gu_chol <- function(x, take.chol = !is.gu_chol(x), ...) {
  out <- if(take.chol) chol(x, ...) else x
  class(out) <- c("gu_chol", class(out)[class(out) != "gu_chol"])
  out
}

#' @rdname utilities
#' @export
is.gu_chol <- function(x) {
  inherits(x, "gu_chol")
}

#' @rdname utilities
#' @export
chol.gu_chol <- function(x, ..., verbose = FALSE) {
  if(verbose) cat("Not taking cholesky again\n")
  x
}


#' @rdname utilities
#' @export
chol_inv <- function(x, ...) {
  chol2inv(chol(x, ...))
}

#' @rdname utilities
#' @export
chol_mvrnorm <- function(n = 1, mu, Sigma, Precision,
                         A = if(missing(Sigma)) backsolve(chol(Precision), diag(1, nrow(Precision))) else t(chol(Sigma))) {
  ## the Cholesky of the inverse t(chol(Sigma)) is the inverse of the Cholesky backsolve(...), but don't get confused:
  ## they might differ in their upper- vs. lower-triangularity.
  ## What's important is that AA' = Sigma
  p <- nrow(A)
  mu <- check_one_or_all(mu, p)
  out <- mu + A %*% matrix(rnorm(p*n), nrow = p)
  if(n == 1) drop(out) else t(out)
}

#' @rdname utilities
#' @export
logit <- function(p) {
  log(p) - log(1 - p)
}

#' @rdname utilities
#' @export
expit <- function(x) {
  1 / (1 + exp(-x))
}



gu_params <- function(...) {
  out <- list(...)
  class(out) <- "gibbs_utils_params"
  out
}

#' @export
print.gibbs_utils_params <- function(x, ...) {
  cat("A 'gibbs_utils_params' object, with:\n")
  for(nm in names(x)) {
    if(is.matrix(x[[nm]])) {
      d <- dim(x[[nm]])
      cat(nm, ": a ", d[1], "x", d[2], " matrix\n", sep = "")
    } else {
      l <- length(x[[nm]])
      if(l == 1 && !is.list(x[[nm]])) {
        cat(nm, ": ", x[[nm]], "\n", sep = "")
      } else {
        cat(nm, ": a length-", l, " vector\n", sep = "")
      }
    }
  }
  invisible(x)
}

