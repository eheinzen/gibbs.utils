% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{utilities}
\alias{utilities}
\alias{gu_chol}
\alias{is.gu_chol}
\alias{chol.gu_chol}
\alias{chol_inv}
\alias{logit}
\alias{expit}
\title{Other Utilities}
\usage{
gu_chol(x, take.chol = !is.gu_chol(x), ...)

is.gu_chol(x)

\method{chol}{gu_chol}(x, ..., verbose = FALSE)

chol_inv(x, ...)

logit(p)

expit(x)
}
\arguments{
\item{x}{A matrix for \code{chol_inv} or a numeric vector for \code{expit}}

\item{take.chol}{Should the cholesky be taken?}

\item{...}{Other arguments (not in use at this time)}

\item{verbose}{Logical, indicating whether to print some stuff.}

\item{p}{A numeric vector}
}
\description{
Other Utilities
}
