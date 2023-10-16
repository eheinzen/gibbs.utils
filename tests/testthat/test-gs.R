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

test_that("diag=TRUE works", {
  X <- cbind(1, c(1, 1, 0, 0))
  Y <- matrix(1:8, 4, 2)
  V <- diag(1, 2)
  Q0 <- diag(1, 4)
  set.seed(20230828)
  conj <- conj_matlm_beta(y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = FALSE)
  set.seed(20230828)
  conjd <- conj_matlm_beta(y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = TRUE)
  expect_equal(conj, conjd)

  set.seed(20230828)
  conj <- gs_matlm_beta(beta = rep(0, 4), y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = FALSE)
  set.seed(20230828)
  conjd <- gs_matlm_beta(beta = rep(0, 4), y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = TRUE)
  expect_equal(conj, conjd)


  ## when there is just one column, the conjugacy should give the same result as the gibbs sampling
  X <- matrix(1, nrow = 4, ncol = 1)
  Y <- matrix(1:8, 4, 2)
  V <- diag(1, 2)
  Q0 <- diag(1, 2)
  set.seed(20230828)
  conj1 <- conj_matlm_beta(y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = FALSE)
  set.seed(20230828)
  conjd1 <- conj_matlm_beta(y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = TRUE)

  set.seed(20230828)
  conj2 <- gs_matlm_beta(beta = rep(0, 2), y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = FALSE)
  set.seed(20230828)
  conjd2 <- gs_matlm_beta(beta = rep(0, 2), y = Y, X = X, V = V, U = NULL, mu0 = NULL, Q0 = Q0, diag = TRUE)
  expect_equal(conj1, conjd1)
  expect_equal(conj2, conjd2)
  expect_equal(conj1, conj2)

})
