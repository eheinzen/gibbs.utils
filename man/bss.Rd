% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bss.R
\name{bss}
\alias{bss}
\alias{B1}
\alias{B2}
\alias{B4}
\title{Calculate Bayesian Smoothing Spline (BSS) for Integers}
\usage{
B1(x)

B2(x)

B4(x)

bss(x, n_basis, first_two = TRUE, from = min(x), to = max(x))
}
\arguments{
\item{x}{for \code{bss}, an integer vector to calculate the basis expansion for.
For \code{B1}, \code{B2}, \code{B4}: a numeric vector between 0 and 1.}

\item{n_basis}{Number of basis terms to use}

\item{first_two}{Should the first two BSS terms (linear and quadratic) be included?}

\item{from, to}{The min and max values that \code{x} can attain.}
}
\description{
Calculate Bayesian Smoothing Spline (BSS) for Integers
}
\details{
Currently this does no linear interpolation and so relies on integer (discrete)
inputs.
}
\references{
Curtis B. Storlie et al. (2015) Calibration of Computational Models
 With Categorical Parameters and Correlated Outputs via Bayesian Smoothing
 Spline ANOVA, Journal of the American Statistical Association, 110:509, 68-82,
 DOI: 10.1080/01621459.2014.979993
}
