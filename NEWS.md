# v1.0.0

- Overhauled the conjugate normal functions, to use `spam::rmvnorm.canonical()` instead of taking two inverses. 
  This method is slower than what we had before when inputs are really small (when the conjugacy is super fast anyway, so I wasn't worried)
  but almost 2x faster when inputs are big. Furthermore, we can get more performance when inputs are `spam` objects already.
  
- `spam` was added as a dependency.

# v0.10.9

- Replaced `crossprod(x, y)` with `t(x) %*% y` in almost all instances, to be compatible with, e.g., the `spam` package for sparse matrices.

- Fixed a typo in the vignette

# v0.10.8

- Added arguments for pre-computed determinants to `dmvnorm()` and `dmatnorm()`.

# v0.10.7

- Applied the same vectorization trick to `dmatnorm()`, for speed gains of over 1000x for large inputs. The default
  for `U=` is also now `NULL`, indicating the identity matrix.

# v0.10.6

- Removed the `use_trace=` argument introduced in v0.9.3. When `x` was short, the trace was
  faster. When `x=` was long, the non-trace method was faster. Now, using a vectorization trick
  with the trace, `dmvnorm()` is always at least as fast (and usually faster) than either of the previous methods.

# v0.10.5

- Also return `V.inv` for conjugate Wishart parameters, `tau` for conjugate normal, and `Q` for some other
  normal conjugate cases where `newQ.inv=` has not already been supplied.
  
- Changed the defaults (but not the default *behavior*) for some normal conjugates.

# v0.10.4

- Changed the default for the `use_trace=` argument introduced in v0.9.3, because in my experience it's faster not to use the matrix trace.

# v0.10.3

- Add the same `mu0 = NULL` default to `gs_diagmatlm_beta()`.

# v0.10.2

- Speed  up (slightly) the conjugacy functions using a default of `mu0 = NULL` to indicate zero.
  This also allows you to give infinite precisions when `mu0 = NULL = 0` (which would otherwise
  give `Inf * 0 = NaN`).

# v0.10.1

- Export `mvqt_binom_approx()` and `mvqt_pois_approx()`.

# v0.10.0

- Added `mh_pois_reg()`.

- Added `proposal="mv quadratic taylor"` to `mh_binom_reg()` to do multivariate Taylor's approximations and MH steps.

# v0.9.3

- Added an option to `dmvnorm()` to opt out of the matrix trace.

# v0.9.2

- Sped up `dmvnorm()` by using the matrix trace instead of a for-loop.

# v0.9.1

- Changed the `n == k == 0` sampling to C++ for `ss_binom_reg()` and `mh_binom_reg()`

- `ss_pois_reg()` now accepts `NA` in `k=`, to indicate that no Poisson draw was observed. 

# v0.9.0

- Renamed `conj_matlm_sigma()` to `conj_matlm_V()`. It, along with `conj_matnorm_V()`, now both use `conj_mvnorm_Q()` in the case that `U = NULL`
  (independence among rows).

# v0.8.2

- Added `truncnorm` as a package dependency.

- Added `gs_diagmatlm_beta()` to Gibbs-sample truncated normal diagonal matrix-normal regression means.

- `cond_mvnorm()` gained optional truncation parameters

- Updated the vignettes to improve the custom LaTeX code that wasn't previously handled by HTML output.

- Improved printing of `"gibbs_utils_params"` objects.

# v0.8.1

- `conj_matnorm_mu()` and `conj_matlm_beta()` gained the `diag=` argument, for 100x speed boost when both `V=` and `Q0=` are diagonal.

- When `params.only=TRUE`, the `mu` vector returned now drops any dimensions in `conj_matnorm_mu()`, `conj_matlm_beta()`, `conj_mvnorm_mu()`,
  `conj_diagmatlm_beta()`, and `conj_lm_beta()`
  
- When `n==1`, `chol_mvrnorm()` now drops dimensions, to match what `MASS::mvrnorm()` does.

- Improved printing of `"gibbs_utils_params"` objects.

# v0.8.0

- `ss_multinom_reg()` now also uses normal draws in the case where `z == 0` or `n == 0`. This is not done multivariately, but rather
  univariately.
  
- `ss_multinom_reg()` now accepts numeric values for the `ref=` argument, to allow inner "columns" to be the reference.

# v0.7.0

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
  
- Updated `mh_binom_reg()` and `ss_binom_reg()` to take normal draws in the multivariate case when only one `n=` is zero instead of
  only doing it when all the `n=`'s are zero.
  
- `mh_binom_reg()` now enforces `k <= n`

- Fixed one critical bug in `mh_binom_reg()` in which the log-likelihood evaluations were using `k=` instead of `n=`.

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
