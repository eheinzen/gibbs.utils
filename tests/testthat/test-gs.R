test_that("Gibbs works for 1-D cases", {
  # this example also appears in test-conjugacy.R, but I have it split out here for consistency
  mu1 <- conj_norm_mu(y = 1, tau = 10, mu0 = 0, tau0 = 1, params.only = TRUE)
  one <- matrix(1, 1, 1)
  ten <- matrix(10, 1, 1)
  mu7 <- conj_matlm_beta(y = one, X = one, V = ten, U = NULL, mu0 = 0, Q0 = one, zero = one, params.only = TRUE)
  mu8 <- gs_matlm_beta(beta = 1.1, y = one, X = one, V = ten, U = NULL, mu0 = 0, Q0 = one, zero = one, params.only = TRUE)

  expect_equal(mu1$mu, mu7$mu)
  expect_equal(mu1$mu, mu8$mu[[1]])
})

test_that("Conjugacy works for 2-D independent cases", {
  mu1 <- c(
    conj_norm_mu(y = 1, tau = 10, mu0 = 0, tau0 = 1, params.only = TRUE)$mu,
    conj_norm_mu(y = 3, tau = 1, mu0 = 2, tau0 = 2, params.only = TRUE)$mu
  )
  one <- matrix(c(1, 3), 1, 2)
  two <- diag(c(1, 2))
  ten <- diag(c(10, 1))
  zero <- c(0, 2)
  X <- matrix(1, 1, 1)
  mu7 <- conj_matlm_beta(y = one, X = cbind(X, X), V = ten, U = NULL, mu0 = zero, Q0 = two, zero = diag(2), params.only = TRUE)
  mu8 <- gs_matlm_beta(beta = c(1.1, 1.65), y = one, X = cbind(X, X), V = ten, U = NULL, mu0 = zero, Q0 = two, zero = diag(2), params.only = TRUE)
  expect_equal(mu1, unlist(mu8$mu))
  expect_equal(mu1, mu7$mu)
})
