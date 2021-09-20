test_that("Conjugacy works for 1-D cases", {
  mu1 <- conj_norm_mu(y = 1, tau = 10, mu0 = 0, tau0 = 1, params.only = TRUE)
  one <- matrix(1, 1, 1)
  ten <- matrix(10, 1, 1)
  mu2 <- conj_mvnorm_mu(y = 1, Q = ten, mu0 = 0, Q0 = one, params.only = TRUE)
  mu3 <- conj_matnorm_mu(y = one, V = ten, U = NULL, mu0 = 0, Q0 = one, params.only = TRUE)
  mu4 <- conj_matlm_beta(y = one, X = one, V = ten, U = NULL, mu0 = 0, Q0 = one, params.only = TRUE)
  mu5 <- conj_lm_beta(y = one, X = one, tau = 10, mu0 = 0, Q0 = one, params.only = TRUE)
  mu6 <- conj_matnorm_mu(y = one, V = one, U = ten, mu0 = 0, Q0 = one, params.only = TRUE)
  mu7 <- conj_diagmatlm_beta(y = one, X = one, V = ten, U = NULL, mu0 = 0, Q0 = one, params.only = TRUE)

  expect_equal(mu1$mu, as.vector(mu2$mu))
  expect_equal(mu1$mu, as.vector(mu3$mu))
  expect_equal(mu1$mu, as.vector(mu4$mu))
  expect_equal(mu1$mu, as.vector(mu5$mu))
  expect_equal(mu1$mu, as.vector(mu6$mu))
  expect_equal(mu1$mu, as.vector(mu7$mu))
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
  mu2 <- conj_mvnorm_mu(y = one, Q = ten, mu0 = zero, Q0 = two, params.only = TRUE)
  mu3 <- conj_matnorm_mu(y = one, V = ten, U = NULL, mu0 = zero, Q0 = two, params.only = TRUE)
  mu4 <- conj_matlm_beta(y = one, X = X, V = ten, U = NULL, mu0 = zero, Q0 = two, params.only = TRUE)
  mu6 <- conj_matnorm_mu(y = t(one), V = X, U = ten, mu0 = zero, Q0 = two, params.only = TRUE)
  mu7 <- conj_diagmatlm_beta(y = one, X = cbind(X, X), V = ten, U = NULL, mu0 = zero, Q0 = two, params.only = TRUE)

  expect_equal(mu1, as.vector(mu2$mu))
  expect_equal(mu1, as.vector(mu3$mu))
  expect_equal(mu1, as.vector(mu4$mu))
  expect_equal(mu1, as.vector(mu6$mu))
  expect_equal(mu1, as.vector(mu7$mu))
})

