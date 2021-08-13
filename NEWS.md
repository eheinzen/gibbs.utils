# v0.3.2

- Defer to `stats::rWishart()` instead of `rwish()`, since the former is faster.

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
