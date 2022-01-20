# v0.6.0.9003

- Sped up the multivariate binomial, multinomial, and Poisson likelihood evaluations (only completing the square once instead of each time),
  so that
  
    - `ss_pois_reg()` is as much as 6.5x faster
    
    - `ss_multinom_reg()` is as much as 1.6x faster
    
    - `ss_binom_reg()` is as much as 4.5x faster
    
    - `mh_binom_reg()` is as much as 2x faster.
    
  As such, there are no longer C++ functions to evaluate the multivariate log-likelihoods (nor the slice sample functions to do the same);
  instead, the conditional means are computed within `cond_mv_mean()`, and the results passed to the univariate likelihoods.
  
- Overhauled `mh_binom_reg()`. It now takes a `proposal=` argument, to determine how proposals are made. The "normal" method
  is the default for backwards-compatibility, but the "quadratic taylor" gives better results and is about the same speed. The C++ internals
  have also been almost completely rewritten.
  
- `mh_binom_reg()` now enforces `k <= n`

# v0.6.0

- `ss_binom_reg()` now performs conjugate normal updates when `n == 0` (or `rowSums(n) == 0` in the multivariate case).

- `ss_binom_reg()` now enforces `k <= n`

- Fixed a bug in `chol_mvrnorm()`.

# v0.5.1

- Small performance improvements for `conj_mvnorm()`.

# v0.5.0

- Added `conj_binom_p()` for the beta conjugacy for a binomial probability.

- Fixed a bug in `ss_multinom_reg()` when dimensions were length 1.

# v0.4.8

- A bug fix for slice sampling multinomial.

# v0.4.7

- Added some arguments for speed and efficiency, and therefore removed some parameters if `params.only = TRUE`.

# v0.4.6

- Updated vignettes.

- Added `conj_diagmatlm_beta()`.

- Added tests for conjugacy.

# v0.4.5

- Added `ref=` argument to `ss_multinom_reg()`.

# v0.4.4

- Sped up multinomial log-likelihood by avoiding exponentiation altogether for structural zeros.

- Changed the acceptable dimensions of `mean=` and `precision=` in `ss_multinom_reg()`

# v0.4.3

- Sped up sampling (by maybe 15%?) by avoiding repeated subsetting of k and n.

- Added `ss_multinom_reg()`.

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
