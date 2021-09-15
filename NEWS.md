# v0.4.3

- Sped up sampling (by maybe 15%?) by avoiding repeated subsetting of k and n.

# v0.4.2

- Sped up binomial log-likelihood (by a non-trivial amount) by avoiding an exponentiation.

# v0.4.1

- Added a pkgdown site (and README.md).

- Added the `params.only=` argument to conjugacy functions.

- Added `impute_conj_mvnorm()`.

# v0.4.0

- Include an option to use (the new function) `chol_mvrnorm()` to draw multivariate normal samples instead of `MASS::mvrnorm()` 
  (which uses eigen decomposition instead of Cholesky).

# v0.3.2

- Defer to `stats::rWishart()` instead of `rwish()`, since the former is faster.

- Added `mh_binom_reg()` (and proceeded to fix one critical bug in 61de4f70).

- Sped up multivariate log-likelihood.

- Fixed one test.

# v0.3.1

- Updated `ss_pois_reg()` and `ss_binom_reg()` to allow for matrix (multivariate) inputs.

- Added tests.

# v0.3.0

- Updated `ss_binom_reg()` to take in and output logit-p instead of p.

# v0.2.0

- Improved matrix normal conjugacies by using a Kronecker product trick.

# v0.1.2

- Fixed a bug in defaults for `conj_mvnorm_mu()`, where the length was accidentally the number of observations instead of variables.

- Improved documentation.

- Improved checking for matrix arguments.

# v0.1.1

- Export `logit()` and `expit()`

- Rename `my_inv()` to `chol_inv()`.

# v0.1.0

Initial version of package
