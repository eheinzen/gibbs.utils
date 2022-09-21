
#' @useDynLib gibbs.utils
#' @importFrom Rcpp sourceCpp
NULL

check_one_or_all <- function(x, len) {
  if(length(x) == 1) rep_len(x, len) else if(length(x) == len) x else stop("'x' must be of length 1 or ", len)
}

vec <- function(...) {
  op <- options(spam.structurebased=FALSE)
  on.exit(options(op))
  as.vector(...)
}

#' Other Utilities
#'
#' @param x A matrix for \code{chol_inv} or a numeric vector for \code{expit}
#' @param p A numeric vector
#' @param ... Other arguments (not in use at this time)
#' @param verbose Logical, indicating whether to print some stuff.
#' @param take.chol Should the cholesky be taken?
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

